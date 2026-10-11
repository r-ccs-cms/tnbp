#include "tnbp/ptns/function/mcmc_sampling.h"
#include "tnbp/ptns/function/mcmc_group_io.h"
#include "tnbp/framework/graphio.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <random>

using S=std::complex<double>;
using T=gqten::tensor<S>;
using Edges=std::vector<std::pair<int,int>>;
using Lines=std::vector<std::vector<int>>;

Lines read_lines(const std::string& name) {
    std::ifstream f(name);if(!f) throw std::runtime_error("cannot read lines");
    Lines lines;std::string row;
    while(std::getline(f,row)) {
        std::istringstream in(row);std::vector<int> line;int site;
        while(in>>site) line.push_back(site);
        if(!line.empty()) lines.push_back(line);
    }
    return lines;
}
void print_plan(const Lines& lines,const tnbp::MCMCGrouping<int>& groups) {
    for(std::size_t l=0;l<lines.size();++l) {
        std::cout<<"line "<<l<<":";for(auto v:lines[l]) std::cout<<' '<<v;std::cout<<'\n';
        for(const auto& g:groups[l]) {std::cout<<"group "<<l<<":";for(auto v:g) std::cout<<' '<<v;std::cout<<'\n';}
    }
    std::cout<<std::flush;
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    std::string stage="arguments";
    try {
        if(argc<6 || argc>8) throw std::runtime_error("usage: heavyhex51 data-dir manual|auto chi samples virtual-dim(1|2) [uniform] [heatbath]");
        bool uniform=false;int method=0;
        for(int i=6;i<argc;++i) {
            const std::string flag=argv[i];
            if(flag=="uniform") uniform=true;
            else if(flag=="heatbath") method=1;
            else throw std::runtime_error("unknown sampling flag");
        }
        const std::string dir=argv[1],mode=argv[2];
        const int chi=std::stoi(argv[3]),count=std::stoi(argv[4]),dim=std::stoi(argv[5]);
        if((mode!="manual" && mode!="auto") || chi<0 || count<1 || (dim!=1 && dim!=2)) throw std::runtime_error("invalid options");
        stage="graph";auto edges=tnbp::read_edges_from_file<int>(dir+"/heavyhex51.edges",true);
        const auto sites=tnbp::GetSiteIndexFromBond(edges);
        if(sites.size()!=51 || edges.size()!=56) throw std::runtime_error("51-site figure fixture count mismatch");
        Lines lines;tnbp::MCMCGrouping<int> groups;
        if(mode=="manual") {
            lines=read_lines(dir+"/horizontal.lines");
            std::ifstream f(dir+"/horizontal.groups");groups=tnbp::read_mcmc_grouping<int>(f,lines);
        } else {
            // This is the only geometry hint in automatic mode.
            tnbp::bmps_setup_lines(edges,std::vector<int>{29,30,31,32,33},lines);
            // Print the candidate even if the public planner later rejects it.
            groups.resize(lines.size());
            for(std::size_t l=0;l<lines.size();++l) {
                for(auto v:lines[l]) groups[l].push_back({v});
                for(int d:{-1,1}) {
                    if((d<0 && l==0)||(d>0 && l+1==lines.size())) continue;
                    auto side=tnbp::bmps_neighbor_groups(edges,lines[l],lines[d<0 ? l-1 : l+1]);
                    for(std::size_t j=0;j<side.size();++j) groups[l][j].insert(groups[l][j].end(),side[j].begin()+1,side[j].end());
                }
            }
        }
        if(rank==0) {std::cout<<"mode="<<mode<<" ranks="<<size<<" D="<<dim<<" chi="<<chi<<" method="<<method<<" uniform="<<uniform<<" samples="<<count<<'\n';print_plan(lines,groups);}
        stage="plan validation";
        if(mode=="manual") tnbp::validate_mcmc_grouping(edges,lines,groups);
        else if(tnbp::make_mcmc_grouping(edges,lines)!=groups) throw std::runtime_error("candidate/public planner differ");
        if(rank==0) std::cout<<"PLAN PASS\n"<<std::flush;
        if(chi==0) {MPI_Finalize();return 0;}
        if(uniform && dim!=1) throw std::runtime_error("uniform control requires D=1");
        stage="state construction";
        tcapi::context_handle_t<T> ctx;tcapi::create_context(ctx);
        std::map<int,int> owners;std::map<int,double> p;
        for(std::size_t i=0;i<sites.size();++i) {owners[sites[i]]=i%size;p[sites[i]]=uniform ? 0.5 : 0.25+0.08*(i%5);}
        // Exact PEPS for product Ry preparations, one RZZ(gamma) on every edge,
        // then Ry(beta) on every site. D=1 is the product-state control (gamma=0).
        const double gamma=dim==2 ? 0.24 : 0.,beta=uniform ? 0. : 0.37;
        const S b[2]={std::sqrt(std::cos(gamma/2)),std::sqrt(S(0,-std::sin(gamma/2)))};
        std::vector<T> W;std::vector<int> ids;
        for(auto site:sites) if(owners.at(site)==rank) {
            const auto degree=tnbp::GetSurroundingBondIndex(site,edges).size();
            tcapi::shape_t<T> shape(degree+1,dim);shape.back()=2;
            auto tensor=tcapi::zeros<T>(ctx,shape);
            const int combinations=dim==2 ? (1<<degree) : 1;
            for(int mask=0;mask<combinations;++mask) for(int s=0;s<2;++s) {
                tcapi::elem_coors_t<T> coordinate;for(std::size_t k=0;k<degree;++k) coordinate.push_back((mask>>k)&1);
                coordinate.push_back(s);S value=0.;
                for(int t=0;t<2;++t) {
                    const double rotation=s==t ? std::cos(beta/2) : (s==1 ? 1. : -1.)*std::sin(beta/2);
                    S term=rotation*std::sqrt(t ? p.at(site) : 1-p.at(site));
                    for(std::size_t k=0;k<degree;++k) term*=b[(mask>>k)&1]*(((mask>>k)&1) && t ? -1. : 1.);
                    value+=term;
                }
                tcapi::set_elem(ctx,tensor,coordinate,value);
            }
            ids.push_back(site);W.push_back(std::move(tensor));
        }
        tnbp::MCMCWorkspace<T,int> work;work.lines=lines;if(mode=="manual") work.manual_grouping=groups;
        std::mt19937 rng(7129+rank);std::vector<int> initial(ids.size(),0);
        const auto start=std::chrono::steady_clock::now();
        stage="boundary initialization";
        tnbp::mcmc_initialize(ctx,edges,W,ids,owners,std::vector<std::vector<int>>{initial},
            10,1.e-9,tcapi::bond_dim_t<T>(chi),1.e-12,1.e-10,work,MPI_COMM_WORLD);
        if(work.site_line_sets!=groups) throw std::runtime_error("stored plan differs");
        if(rank==0) std::cout<<"INITIALIZATION PASS\n"<<std::flush;
        const int burn=16;
        std::vector<std::vector<std::vector<int>>> samples(burn+count,std::vector<std::vector<int>>(1));
        samples[0][0]=initial;stage="sampling";
        tnbp::mcmc_sampling(ctx,edges,W,ids,owners,samples,method,rng,
            10,1.e-9,tcapi::bond_dim_t<T>(chi),1.e-12,1.e-10,work,MPI_COMM_WORLD);
        std::vector<int> local((burn+count)*51,0),global(local.size());
        for(int n=0;n<burn+count;++n) {
            if(samples[n][0].size()!=ids.size()) throw std::runtime_error("wrong sample length");
            for(std::size_t j=0;j<ids.size();++j) {
                const auto s=samples[n][0][j];if(s!=0 && s!=1) throw std::runtime_error("nonbinary sample");
                const auto k=std::distance(sites.begin(),std::find(sites.begin(),sites.end(),ids[j]));local[n*51+k]=s;
            }
        }
        MPI_Allreduce(local.data(),global.data(),int(local.size()),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        if(rank==0) {
            double max_error=0,mean_observed=0,mean_expected=0;
            for(std::size_t j=0;j<sites.size();++j) {
                S phase=1.;const auto site=sites[j];
                for(auto [a,c]:edges) if(a==site || c==site) {
                    const double q=p.at(a==site ? c : a);
                    phase*=(1-q)*std::exp(S(0,-gamma))+q*std::exp(S(0,gamma));
                }
                const double z=(1-2*p.at(site))*std::cos(beta)-2*std::sqrt(p.at(site)*(1-p.at(site)))*std::sin(beta)*phase.real();
                const double expected=(1-z)/2;double observed=0.;
                for(int n=burn;n<burn+count;++n) observed+=global[n*51+j]/double(count);
                max_error=std::max(max_error,std::abs(expected-observed));mean_observed+=observed;mean_expected+=expected;
                std::cout<<"marginal "<<site<<" expected="<<expected<<" observed="<<observed<<'\n';
            }
            for(int n=burn;n<burn+count;++n) {std::cout<<"sample "<<n-burn<<": ";for(int j=0;j<51;++j) std::cout<<global[n*51+j];std::cout<<'\n';}
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            std::cout<<"SAMPLING COMPLETED mean_ones="<<mean_observed<<" exact_mean_ones="<<mean_expected
                     <<" max_marginal_error="<<max_error<<" seconds="<<seconds<<'\n';
        }
        tcapi::destroy_context(ctx);
    } catch(const std::exception& e) {
        std::cerr<<"FAILED rank="<<rank<<" stage="<<stage<<": "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    MPI_Finalize();return 0;
}
