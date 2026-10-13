#include "../common/device.h"
#include "tnbp/ptns/function/mcmc_sampling.h"
#include "tnbp/ptns/function/mcmc_group_io.h"
#include <array>
#include <iostream>
#include <random>
#include <stdexcept>

using Scalar = std::complex<double>;
#ifdef TNBP_TEST_CUDA
using Tensor = tcapi::cuda::Tensor<Scalar>;
#else
using Tensor = gqten::tensor<Scalar>;
#endif
using Context = tcapi::context_handle_t<Tensor>;
using Edges = std::vector<std::pair<int,int>>;

// Independent dense reference: enumerate virtual indices, without TNBP contraction.
#include "../common/dense_fixture.h"
Scalar exact(const Edges& edges, int state) {
    const auto nsites=tnbp::GetSiteIndexFromBond(edges).size();
    Scalar amplitude = 0.;
    for (int virtuals=0; virtuals<(1<<edges.size()); ++virtuals) {
        Scalar product = 1.;
        for (int site=0; site<int(nsites); ++site) {
            std::vector<int> bonds;
            for (std::size_t e=0; e<edges.size(); ++e)
                if(edges[e].first==site || edges[e].second==site)
                    bonds.push_back((virtuals>>e)&1);
            product *= entry(site,(state>>site)&1,bonds);
        }
        amplitude += product;
    }
    return amplitude;
}
void check_bmps(Context& ctx,const Edges& edges,const std::vector<Tensor>& attached,
                const std::vector<int>& ids,const std::map<int,int>& owners,
                const std::vector<std::vector<int>>& lines,
                const tnbp::MCMCGrouping<int>* grouping,Scalar expected,MPI_Comm comm) {
    auto graph_only=edges,tensor_edges=edges;
    auto local_ids=ids;auto ranks=owners;
    std::vector<Tensor> tensors;
    for(const auto& t:attached) tensors.push_back(tcapi::copy(ctx,t));
    const auto* groups=grouping ? &(*grouping)[1] : nullptr;
    tnbp::bmps_increment_line_for_edges(lines[0],lines[1],graph_only,groups);
    tnbp::bmps_increment_line(ctx,lines[0],lines[1],tensor_edges,tensors,local_ids,ranks,comm,groups);
    if(graph_only!=tensor_edges) throw std::runtime_error("BMPS tensor/graph plans disagree");
    for(std::size_t j=0;j<lines[1].size();++j) {
        if(ranks.at(lines[1][j])!=owners.at(lines[1][j])) throw std::runtime_error("BMPS shifted owner");
        for(std::size_t k=j+2;k<lines[1].size();++k)
            if(tnbp::find_edge(tensor_edges,lines[1][j],lines[1][k])!=tensor_edges.end())
                throw std::runtime_error("BMPS absorption retained a long bond");
    }
    for(int site:lines[0]) if(ranks.count(site)) throw std::runtime_error("BMPS retained absorbed site");
    const auto schedule=tensor_edges;
    tnbp::graph_tensor_network_contraction(ctx,tensor_edges,tensors,local_ids,ranks,
                                          schedule,std::map<int,int>{},std::map<int,int>{},comm);
    if(!tensor_edges.empty() || ranks.size()!=1) throw std::runtime_error("BMPS final graph");
    int rank;MPI_Comm_rank(comm,&rank);
    if(rank==ranks.begin()->second) {
        if(tensors.size()!=1 || tcapi::order(ctx,tensors[0])!=0 ||
           std::abs(tcapi::get_elem(ctx,tensors[0],{})-expected)>1.e-10*(1.+std::abs(expected)))
            throw std::runtime_error("BMPS absorption changes exact amplitude");
    }
}
std::vector<int> run(Context& ctx, const Edges& edges, const std::vector<std::vector<int>>& lines,
         const std::string& name, MPI_Comm comm, const tnbp::MCMCGrouping<int>* grouping=nullptr, int fixed_mask=0, bool force_workspace=false, int method=0) {
    int rank, size; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
    const auto nsites=tnbp::GetSiteIndexFromBond(edges).size();
    const int nstates=1<<nsites;
    const bool workspace_mode=grouping || force_workspace;
    std::vector<Tensor> W;
    std::vector<int> ids;
    std::map<int,int> owners;
    for(int site=0; site<int(nsites); ++site) {
        owners[site]=site%size;
        if(site%size!=rank) continue;
        ids.push_back(site);
        const auto bonds=tnbp::GetSurroundingBondIndex(site,edges);
        tcapi::shape_t<Tensor> shape(bonds.size()+1,2);
        const int physical_dim=(fixed_mask & (1<<site)) ? 1 : 2;
        shape.back()=physical_dim;
        auto t=tcapi::zeros<Tensor>(ctx,shape);
        for(int b=0; b<(1<<bonds.size()); ++b) for(int s=0;s<physical_dim;++s) {
            tcapi::elem_coors_t<Tensor> coors;
            std::vector<int> values;
            for(std::size_t k=0;k<bonds.size();++k) {
                coors.push_back((b>>k)&1); values.push_back((b>>k)&1);
            }
            coors.push_back(s);
            tcapi::set_elem(ctx,t,coors,entry(site,s,values));
        }
        W.push_back(std::move(t));
    }
    std::vector<double> probability(nstates);
    double norm=0.;
    for(int state=0;state<nstates;++state) {
        if(state & fixed_mask) continue;
        const Scalar expected=exact(edges,state);
        probability[state]=std::norm(expected); norm+=probability[state];
        std::vector<int> local_states;
        for(int site:ids) local_states.push_back((state>>site)&1);
        std::vector<Tensor> attached;
        tnbp::mcmc_init_attachment(ctx,edges,W,ids,owners,local_states,lines,attached);
        check_bmps(ctx,edges,attached,ids,owners,lines,grouping,expected,comm);
        auto remaining=edges; auto local_ids=ids; auto ranks=owners;
        std::map<int,int> priority{{0,4},{1,3},{2,2},{3,1}}, labels;
        tnbp::graph_tensor_network_contraction(ctx,remaining,attached,local_ids,ranks,
                                             edges,priority,labels,comm);
        if(!remaining.empty() || ranks.size()!=1 || ranks.count(0)!=1)
            throw std::runtime_error(name+": contraction metadata");
        if(rank==owners.at(0)) {
            if(attached.size()!=1 || tcapi::order(ctx,attached[0])!=0 ||
               std::abs(tcapi::get_elem(ctx,attached[0],{})-expected)>1.e-10*(1.+std::abs(expected)))
                throw std::runtime_error(name+": exact amplitude mismatch");
        }
    }
    if(rank==0) std::cout<<name<<": exact contraction PASS; starting boundary initialization"<<std::endl;
    for(double& p:probability) p/=norm;
    std::vector<std::vector<Tensor>> V(1);
    std::vector<Edges> line_edges;
    std::vector<std::vector<std::vector<Tensor>>> upper(1), lower(1), eu(1), el(1);
    std::vector<std::vector<int>> iu,il,ie_u,ie_l;
    std::vector<std::vector<std::vector<int>>> sets;
    std::vector<int> initial(ids.size(),0);
    tnbp::MCMCWorkspace<Tensor,int> workspace;
    if(workspace_mode) {
        workspace.lines=lines;
        if(grouping) workspace.manual_grouping=*grouping;
        tnbp::mcmc_initialize(ctx,edges,W,ids,owners,std::vector<std::vector<int>>{initial},
            30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,workspace,comm);
    } else tnbp::mcmc_initialize(ctx,edges,W,ids,owners,initial,lines,V[0],line_edges,
        upper[0],iu,eu[0],ie_u,lower[0],il,el[0],ie_l,sets,
        30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,comm,grouping);
    // Both directions return only proper prefixes/suffixes, with original site
    // labels. Seed output buffers to also check clearing for a one-site line.
    const auto& graphs=workspace_mode ? workspace.line_edges : line_edges;
    const auto& blocks=workspace_mode ? workspace.site_line_sets : sets;
    if(!grouping && blocks!=tnbp::make_mcmc_grouping(edges,lines))
        throw std::runtime_error("initialization did not retain automatic plan");
    const auto& values=workspace_mode ? workspace.V[0] : V[0];
    const auto& up=workspace_mode ? workspace.T_upper[0] : upper[0];
    const auto& down=workspace_mode ? workspace.T_lower[0] : lower[0];
    const auto& up_ids=workspace_mode ? workspace.site_idx_upper : iu;
    const auto& down_ids=workspace_mode ? workspace.site_idx_lower : il;
    const auto& down_messages=workspace_mode ? workspace.E_lower[0] : el[0];
    const auto& down_message_ids=workspace_mode ? workspace.edge_idx_lower : ie_l;
    if(!down.back().empty() || !down_ids.back().empty() ||
       !down_messages.back().empty() || !down_message_ids.back().empty())
        throw std::runtime_error("unused uppermost lower boundary was initialized");
    const std::vector<int> empty_ids;
    const std::vector<Tensor> empty_tensors;
    for(std::size_t l=0;l<lines.size();++l) {
        std::vector<int> target_ids;std::vector<Tensor> target;
        for(std::size_t j=0;j<ids.size();++j)
            if(std::find(lines[l].begin(),lines[l].end(),ids[j])!=lines[l].end()) {
                target_ids.push_back(ids[j]);target.push_back(tcapi::copy(ctx,values[j]));
            }
        for(int direction:{0,1}) {
            std::vector<Tensor> env;
            env.push_back(tcapi::zeros<Tensor>(ctx,tcapi::shape_t<Tensor>{}));
            std::vector<int> env_ids{-999},expected;
            for(std::size_t step=0;step+1<lines[l].size();++step) {
                const auto site=lines[l][direction==0 ? step : lines[l].size()-1-step];
                if(owners.at(site)==rank) expected.push_back(site);
            }
            tnbp::mcmc_init_environment(ctx,owners,graphs[l],
                l+1<lines.size() ? lines[l+1] : empty_ids,
                l+1<lines.size() ? up[l+1] : empty_tensors,
                l+1<lines.size() ? up_ids[l+1] : empty_ids,
                lines[l],target,target_ids,
                l ? lines[l-1] : empty_ids,l ? down[l-1] : empty_tensors,
                l ? down_ids[l-1] : empty_ids,env,env_ids,direction,blocks[l],comm);
            if(env_ids!=expected || env.size()!=expected.size())
                throw std::runtime_error("environment prefix/suffix count or labels");
        }
    }
    if(rank==0) std::cout<<name<<": boundary initialization PASS; starting sampling"<<std::endl;
    const int burn=fixed_mask==nstates-1 ? 4 : 128, count=fixed_mask==nstates-1 ? 64 : 4096;
    std::vector<std::vector<std::vector<int>>> samples(burn+count);
    for(auto& sample:samples) sample.resize(1);
    samples[0][0]=initial;
    std::mt19937 rng(1729+rank);
    if(workspace_mode) {
        tnbp::mcmc_sampling(ctx,edges,W,ids,owners,samples,method,rng,
            30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,workspace,comm);
    } else tnbp::mcmc_sampling(ctx,edges,W,ids,owners,lines,samples,line_edges,V,
        upper,iu,eu,ie_u,lower,il,el,ie_l,sets,method,rng,
        30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,comm,grouping);
    std::vector<int> local(samples.size()),global(samples.size());
    for(std::size_t n=0;n<samples.size();++n) {
        if(samples[n][0].size()!=ids.size()) throw std::runtime_error("sample size");
        for(std::size_t j=0;j<ids.size();++j) {
            const int s=samples[n][0][j];
            if(s<0||s>1 || ((fixed_mask & (1<<ids[j])) && s!=0))
                throw std::runtime_error("sample state (including fixed sites)");
            local[n]|=s<<ids[j];
        }
    }
    MPI_Allreduce(local.data(),global.data(),int(local.size()),MPI_INT,MPI_SUM,comm);
    std::vector<int> histogram(nstates);
    for(int n=burn;n<burn+count;++n) ++histogram.at(global[n]);
    double max_error=0.,total_variation=0.;
    for(int s=0;s<nstates;++s) {
        const auto error=std::abs(double(histogram[s])/count-probability[s]);
        max_error=std::max(max_error,error);total_variation+=0.5*error;
    }
    if(max_error>0.04 || total_variation>0.10) throw std::runtime_error(name+": sampled distribution mismatch");
    if(rank==0) std::cout<<name<<": allowed-state exact amplitudes PASS; "<<count<<" samples max error="<<max_error<<" TV="<<total_variation<<" PASS\n";
    return global;
}
void plan_checks() {
    const Edges edges{{0,1},{0,2},{0,3},{1,3},{2,3}};
    const std::vector<std::vector<int>> lines{{0,1},{2,3}};
    std::istringstream input("# paired groups\n0 0 2\n0 1 3\n1 2 0\n1 3 1\n");
    const auto groups=tnbp::read_mcmc_grouping<int>(input,lines);
    tnbp::validate_mcmc_grouping(edges,lines,groups);
    const auto schedule=tnbp::bmps_absorption_schedule(lines[0],lines[1],edges,&groups[1]);
    if(schedule!=Edges{{2,0},{3,1}}) throw std::runtime_error("manual absorption ignored");
    const Edges diagonal12{{0,1},{0,2},{1,2},{1,3},{2,3}};
    const auto automatic=tnbp::bmps_absorption_schedule(lines[0],lines[1],diagonal12);
    const auto manual=tnbp::bmps_absorption_schedule(lines[0],lines[1],diagonal12,&groups[1]);
    if(automatic!=manual || manual!=Edges{{2,0},{3,1}})
        throw std::runtime_error("BMPS automatic/manual policy selection");
    for(const auto& graph:{edges,diagonal12}) {
        const auto automatic_groups=tnbp::make_mcmc_grouping(graph,lines);
        if(automatic_groups!=groups) throw std::runtime_error("nearest unassigned diagonal selection");
        auto reversed=graph;std::reverse(reversed.begin(),reversed.end());
        for(auto& e:reversed) std::swap(e.first,e.second);
        if(tnbp::make_mcmc_grouping(reversed,lines)!=automatic_groups)
            throw std::runtime_error("grouping depends on edge order");
    }
    const std::vector<std::vector<int>> sparse_lines{{0,1,2,3,4,5},{6,7,8,9,10,11}};
    Edges sparse;
    for(int k=0;k<5;++k) {sparse.emplace_back(k,k+1);sparse.emplace_back(k+6,k+7);}
    for(int k:{0,2,4}) sparse.emplace_back(k,k+6);
    const auto sparse_groups=tnbp::make_mcmc_grouping(sparse,sparse_lines);
    if(sparse_groups[1]!=std::vector<std::vector<int>>{{6,0,1},{7},{8,2,3},{9},{10,4,5},{11}})
        throw std::runtime_error("sparse rung coverage/distance priority");
    // Relabeling must not replace the specified line order as tie breaker.
    auto relabeled=sparse;auto relabeled_lines=sparse_lines;
    const std::vector<int> labels{91,4,70,2,35,8,19,6,54,1,43,7};
    for(auto& e:relabeled) {e.first=labels[e.first];e.second=labels[e.second];}
    for(auto& line:relabeled_lines) for(auto& site:line) site=labels[site];
    auto expected_groups=sparse_groups;
    for(auto& window:expected_groups) for(auto& group:window) for(auto& site:group) site=labels[site];
    if(tnbp::make_mcmc_grouping(relabeled,relabeled_lines)!=expected_groups)
        throw std::runtime_error("automatic tie break uses labels instead of line order");
    auto rejects=[](auto operation) {
        bool rejected=false;
        try { operation(); } catch(const std::invalid_argument&) { rejected=true; }
        if(!rejected) throw std::runtime_error("invalid manual plan accepted");
    };
    rejects([&]{tnbp::make_mcmc_grouping(Edges{{0,1},{2,3}},lines);});
    const Edges unsupported{{0,1},{0,2},{1,4},{2,3},{2,5},{4,5},{3,5}};
    const std::vector<std::vector<int>> unsupported_lines{{0,1},{2,3,4},{5}};
    rejects([&]{tnbp::make_mcmc_grouping(unsupported,unsupported_lines);});
    tnbp::validate_mcmc_grouping(unsupported,unsupported_lines,
        tnbp::MCMCGrouping<int>{{{0,2,3},{1,4}},{{2,0},{3},{4,1,5}},{{5,2,4,3}}});
    auto unchanged=edges;
    rejects([&]{tnbp::bmps_increment_line_for_edges(lines[0],std::vector<int>{},unchanged);});
    if(unchanged!=edges) throw std::runtime_error("BMPS invalid plan mutated graph");
    rejects([&]{tnbp::bmps_absorption_schedule(std::vector<int>{99},lines[1],edges);});
    auto bad=groups;bad[0][1].push_back(2);
    rejects([&]{tnbp::validate_mcmc_grouping(edges,lines,bad);});
    bad=groups;bad[0][1].pop_back();
    rejects([&]{tnbp::validate_mcmc_grouping(edges,lines,bad);});
    bad=groups;bad[0][1][1]=99;
    rejects([&]{tnbp::validate_mcmc_grouping(edges,lines,bad);});
    bad=groups;bad[0][0]={0,1,2};
    rejects([&]{tnbp::validate_mcmc_grouping(edges,lines,bad);});
    rejects([&]{std::istringstream x("0 0 2x\n");tnbp::read_mcmc_grouping<int>(x,lines);});
    rejects([&]{std::istringstream x("0 0 2\n0 0 3\n");tnbp::read_mcmc_grouping<int>(x,lines);});
    rejects([&]{std::istringstream x("0 0 2\n");tnbp::read_mcmc_grouping<int>(x,lines);});
    auto disconnected=edges;
    disconnected.erase(disconnected.begin()+1); // remove 0--2
    rejects([&]{tnbp::validate_mcmc_grouping(disconnected,lines,groups);});
    // Explicit relabeling must work without contraction and must preserve edge order.
    auto renamed=edges;std::vector<int> ids{0,1,2,3};std::map<int,int> owners{{0,0},{1,0},{2,0},{3,0}};
    tnbp::graph_rename_site(renamed,ids,owners,1,5);
    if(ids!=std::vector<int>{0,5,2,3} || owners.count(1) || !owners.count(5) ||
       renamed!=Edges{{0,5},{0,2},{0,3},{3,5},{2,3}})
        throw std::runtime_error("singleton relabel corrupts metadata");
    rejects([&]{tnbp::graph_rename_site(renamed,ids,owners,0,5);});
}
void invalid_line_checks(Context& ctx,MPI_Comm comm) {
    const Edges graph{{0,1},{1,2},{2,3}};
    const std::vector<std::vector<std::vector<int>>> invalid{
        {},{{0,1,2,3}},{{0,1},{},{2,3}},{{0,1},{1,2,3}},
        {{0,1},{2}},{{0,1},{2,99}},{{0,2},{1},{3}}};
    const std::vector<Tensor> tensors;
    const std::vector<int> ids,states;
    const std::map<int,int> owners;
    auto rejects=[](auto operation) {
        try {operation();} catch(const std::invalid_argument&) {return;}
        throw std::runtime_error("invalid sampling lines were not rejected");
    };
    for(const auto& lines:invalid) {
        rejects([&]{tnbp::validate_mcmc_lines(graph,lines);});
        // Empty inputs deliberately prove validation precedes all tensor access.
        tnbp::MCMCWorkspace<Tensor,int> w;w.lines=lines;
        w.V.resize(1);w.T_upper.resize(1);w.T_lower.resize(1);
        w.E_upper.resize(1);w.E_lower.resize(1);
        for(bool manual:{false,true}) {
            const tnbp::MCMCGrouping<int> groups;
            rejects([&]{tnbp::mcmc_initialize(ctx,graph,tensors,ids,owners,states,lines,
                w.V[0],w.line_edges,w.T_upper[0],w.site_idx_upper,w.E_upper[0],w.edge_idx_upper,
                w.T_lower[0],w.site_idx_lower,w.E_lower[0],w.edge_idx_lower,w.site_line_sets,
                30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,comm,
                manual ? &groups : nullptr);});
            std::vector<std::vector<std::vector<int>>> samples;
            std::mt19937 rng(1);
            rejects([&]{tnbp::mcmc_sampling(ctx,graph,tensors,ids,owners,lines,samples,
                w.line_edges,w.V,w.T_upper,w.site_idx_upper,w.E_upper,w.edge_idx_upper,
                w.T_lower,w.site_idx_lower,w.E_lower,w.edge_idx_lower,w.site_line_sets,
                0,rng,30,1.e-12,tcapi::bond_dim_t<Tensor>(16),1.e-14,1.e-14,comm,
                manual ? &groups : nullptr);});
        }
        rejects([&]{tnbp::mcmc_initialize(ctx,graph,tensors,ids,owners,
            std::vector<std::vector<int>>{},30,1.e-12,tcapi::bond_dim_t<Tensor>(16),
            1.e-14,1.e-14,w,comm);});
        if(w.V.size()!=1 || !w.line_edges.empty())
            throw std::runtime_error("invalid lines mutated workspace");
    }
    int rank;MPI_Comm_rank(comm,&rank);
    if(rank==0) std::cout<<"invalid line checks PASS"<<std::endl;
}
// A three-state product distribution checks the pairwise (not full Gibbs)
// interpretation, including escape from an initially zero-weight state.
void heat_bath_product_check(Context& ctx, MPI_Comm comm, bool uniform_binary) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const Edges edges{{0,1}};
    double p[2][3]={{0.,.3,.7},{.2,.3,.5}};
    const int dimension=uniform_binary ? 2 : 3;
    if(uniform_binary) for(auto& row:p) {row[0]=row[1]=.5;row[2]=0.;}
    std::map<int,int> owners{{0,0},{1,1%size}};
    std::vector<int> ids;std::vector<Tensor> W;
    for(int site=0;site<2;++site) if(owners.at(site)==rank) {
        ids.push_back(site);
        auto t=tcapi::zeros<Tensor>(ctx,{1,dimension});
        for(int state=0;state<dimension;++state) tcapi::set_elem(ctx,t,{0,state},Scalar(std::sqrt(p[site][state])));
        W.push_back(std::move(t));
    }
    tnbp::MCMCWorkspace<Tensor,int> work;work.lines={{0},{1}};
    std::vector<int> initial(ids.size(),0);
    tnbp::mcmc_initialize(ctx,edges,W,ids,owners,std::vector<std::vector<int>>{initial},
        10,1.e-12,tcapi::bond_dim_t<Tensor>(1),1.e-14,1.e-14,work,comm);
    const int burn=128,count=4096;
    std::vector<std::vector<std::vector<int>>> samples(burn+count,std::vector<std::vector<int>>(1));
    samples[0][0]=initial;std::mt19937 rng(3129+rank);
    tnbp::mcmc_sampling(ctx,edges,W,ids,owners,samples,1,rng,
        10,1.e-12,tcapi::bond_dim_t<Tensor>(1),1.e-14,1.e-14,work,comm);
    std::vector<int> local(count),global(count);
    for(int n=0;n<count;++n) for(std::size_t k=0;k<ids.size();++k) {
        int v=samples[n+burn][0][k];
        if(v<0 || v>=dimension) throw std::runtime_error("product sample out of range");
        local[n]+=v*(ids[k]==0 ? 1 : 3);
    }
    MPI_Allreduce(local.data(),global.data(),count,MPI_INT,MPI_SUM,comm);
    std::array<int,9> histogram{};for(int v:global) ++histogram.at(v);
    for(int b=0;b<3;++b) for(int a=0;a<3;++a)
        if(std::abs(double(histogram[a+3*b])/count-p[0][a]*p[1][b])>.04)
            throw std::runtime_error("ternary heat bath distribution mismatch");
    if(rank==0) std::cout<<(uniform_binary ? "uniform binary" : "ternary zero-weight escape")<<": heat bath distribution PASS\n";
}
int main(int argc,char**argv) {
    MPI_Init(&argc,&argv);
    initialize_tensor_device(MPI_COMM_WORLD);
    try {
        plan_checks();
        Context ctx; tcapi::create_context(ctx);
        invalid_line_checks(ctx,MPI_COMM_WORLD);
        if(argc>1) {
            if(std::string(argv[1])=="heat-bath") {
                for(double scale : {1.e-200,1.,1.e200}) {
                    if(tnbp::mcmc_pair_metropolis_probability(scale,2*scale)!=1. ||
                       tnbp::mcmc_pair_metropolis_probability(2*scale,scale)!=.25 ||
                       tnbp::mcmc_pair_metropolis_probability(scale,scale)!=1.)
                        throw std::runtime_error("Metropolis probability scaling");
                }
                if(tnbp::mcmc_pair_metropolis_probability(0.,1.)!=1. ||
                   tnbp::mcmc_pair_metropolis_probability(1.,0.)!=0. ||
                   tnbp::mcmc_pair_metropolis_probability(0.,0.)!=0.)
                    throw std::runtime_error("Metropolis zero weights");
                for(double scale:{1.e-250,1.,1.e250}) {
                    if(std::abs(tnbp::mcmc_pair_heat_bath_probability(scale,2*scale)-.8)>1.e-14)
                        throw std::runtime_error("heat bath probability scaling");
                }
                if(tnbp::mcmc_pair_heat_bath_probability(0.,1.)!=1. ||
                   tnbp::mcmc_pair_heat_bath_probability(1.,0.)!=0. ||
                   tnbp::mcmc_pair_heat_bath_probability(0.,0.)!=0. ||
                   tnbp::mcmc_pair_heat_bath_probability(1.,1.)!=.5)
                    throw std::runtime_error("heat bath zero/equal weights");
                for(int method:{-1,2}) {
                    bool caught=false;
                    try {tnbp::validate_mcmc_sampling_method(method);}
                    catch(const std::invalid_argument&) {caught=true;}
                    if(!caught) throw std::runtime_error("invalid method accepted");
                }
                heat_bath_product_check(ctx,MPI_COMM_WORLD,false);
                heat_bath_product_check(ctx,MPI_COMM_WORLD,true);
                const Edges graph{{0,1},{1,2},{2,3}};
                const std::vector<std::vector<int>> lines{{0},{1,2},{3}};
                const auto plan=tnbp::make_mcmc_grouping(graph,lines);
                auto a=run(ctx,graph,lines,"heat bath automatic",MPI_COMM_WORLD,nullptr,0,false,1);
                auto b=run(ctx,graph,lines,"heat bath manual",MPI_COMM_WORLD,&plan,0,false,1);
                if(a!=b) throw std::runtime_error("heat bath automatic/manual differ");
                run(ctx,graph,lines,"heat bath fixed",MPI_COMM_WORLD,&plan,15,false,1);
                MPI_Finalize();return 0;
            }
            if(std::string(argv[1])=="automatic-groups") {
                const std::vector<std::vector<int>> lines{{0,1,2},{3}};
                for(bool skip:{false,true}) {
                    Edges graph{{0,1},{1,2},{2,3}};
                    if(skip) graph.emplace_back(0,2);
                    const auto plan=tnbp::make_mcmc_grouping(graph,lines);
                    const auto automatic_samples=run(ctx,graph,lines,skip ? "auto skip" : "auto chain",MPI_COMM_WORLD);
                    const auto manual_samples=run(ctx,graph,lines,skip ? "same manual skip" : "same manual chain",MPI_COMM_WORLD,&plan);
                    if(automatic_samples!=manual_samples) throw std::runtime_error("automatic/manual sampling paths differ");
                }
                Edges sparse;
                for(int k=0;k<5;++k) {sparse.emplace_back(k,k+1);sparse.emplace_back(k+6,k+7);}
                for(int k:{0,2,4}) sparse.emplace_back(k,k+6);
                // All physical sites fixed, virtual dimension 2: exercises the
                // 12-site boundary/environment numerics without 4096 amplitudes.
                run(ctx,sparse,{{0,1,2,3,4,5},{6,7,8,9,10,11}},
                    "automatic sparse workspace",MPI_COMM_WORLD,nullptr,4095,true);
                for(const std::string diagonal:{"diagonal03","diagonal12"}) {
                    Edges graph{{0,1},{0,2},{1,3},{2,3}};
                    graph.push_back(diagonal=="diagonal03" ? std::make_pair(0,3) : std::make_pair(1,2));
                    run(ctx,graph,{{0,1},{2,3}},"automatic "+diagonal,MPI_COMM_WORLD);
                }
            } else if(std::string(argv[1])=="shift-absorption") {
                const std::vector<std::vector<int>> lines{{0,1},{2,3,4},{5}};
                for(bool adjacent:{true,false}) {
                    tnbp::MCMCGrouping<int> groups{{{0,2,3},{1,4}},{{2,0},{3},{4,1,5}},{{5,2,4,3}}};
                    Edges edges{{0,1},{0,2},{1,4},{2,3},{2,5},{4,5}};
                    if(adjacent) edges.push_back({3,4});
                    else {edges.push_back({3,5});std::reverse(edges.begin(),edges.end());for(auto& e:edges) std::swap(e.first,e.second);}
                    run(ctx,edges,lines,adjacent ? "shift adjacent" : "shift absent/reversed",MPI_COMM_WORLD,&groups);
                }
            } else if(std::string(argv[1])=="fixed-sites") {
                tnbp::MCMCGrouping<int> groups{{{0},{1},{2,3}},{{3,2,1,0}}};
                for(int mask:{15,2,5})
                    run(ctx,{{0,1},{1,2},{2,3}},{{0,1,2},{3}},
                        "fixed mask "+std::to_string(mask),MPI_COMM_WORLD,&groups,mask);
            } else if(std::string(argv[1])=="singleton-ends") {
                const Edges edges{{0,1},{2,3},{0,2},{1,3},{0,3}};
                const std::vector<std::vector<int>> lines{{1},{0,3},{2}};
                run(ctx,edges,lines,"singleton ends automatic",MPI_COMM_WORLD);
                tnbp::MCMCGrouping<int> groups{{{1,0,3}},{{0,1,2},{3}},{{2,0,3}}};
                run(ctx,edges,lines,"singleton ends manual",MPI_COMM_WORLD,&groups);
                tnbp::MCMCGrouping<int> three{{{0},{1},{2,3}},{{3,2,1,0}}};
                run(ctx,{{0,1},{1,2},{2,3}},{{0,1,2},{3}},
                    "manual three-site line",MPI_COMM_WORLD,&three);

            } else if(std::string(argv[1])=="manual") {
                tnbp::MCMCGrouping<int> groups{{{0,2},{1,3}},{{2,0},{3,1}}};
                run(ctx,{{0,1},{0,2},{0,3},{1,3},{2,3}},{{0,1},{2,3}},"manual diagonal03",MPI_COMM_WORLD,&groups);
                run(ctx,{{0,1},{0,2},{1,2},{1,3},{2,3}},{{0,1},{2,3}},"manual diagonal12",MPI_COMM_WORLD,&groups);
                tnbp::MCMCGrouping<int> singleton{{{0,2,3},{1}},{{2,0},{3,1}}};
                run(ctx,{{0,1},{0,2},{0,3},{1,3},{2,3}},{{0,1},{2,3}},"manual singleton",MPI_COMM_WORLD,&singleton);
                tnbp::MCMCGrouping<int> chain{{{0,1,2}},{{1,0},{2,3}},{{3,2,1}}};
                run(ctx,{{0,1},{1,2},{2,3}},{{0},{1,2},{3}},"manual three lines",MPI_COMM_WORLD,&chain);
            } else {
            if(std::string(argv[1])!="diagonal03" && std::string(argv[1])!="diagonal12")
                throw std::runtime_error("expected diagonal03 or diagonal12");
            Edges edges{{0,1},{0,2},{1,3},{2,3}};
            edges.push_back(std::string(argv[1])=="diagonal03" ? std::make_pair(0,3) : std::make_pair(1,2));
            std::sort(edges.begin(),edges.end());
            run(ctx,edges,{{0,1},{2,3}},argv[1],MPI_COMM_WORLD);
            }
        } else {
        run(ctx,{{0,1},{1,2},{2,3}},{{0},{1},{2},{3}},"chain",MPI_COMM_WORLD);
        run(ctx,{{0,1},{0,2},{1,3},{2,3}},{{0,1},{2,3}},"square",MPI_COMM_WORLD);
        }
        tcapi::destroy_context(ctx);
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();
}
