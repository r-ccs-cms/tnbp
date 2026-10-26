/// This file is a part of r-ccs-cms/tnbp
#ifndef TNBP_PTNS_FUNCTION_MCMC_EDGE_HEAT_BATH_H
#define TNBP_PTNS_FUNCTION_MCMC_EDGE_HEAT_BATH_H
#include "tnbp/ptns/function/mcmc_heat_bath.h"
#include "tnbp/ptns/function/mcmc_edge_plan.h"
namespace tnbp {
// A retained cut has one axis per original edge, in original edge order.
// Inside a region, tensors stay in native graph representation across steps.
template<class TenT,class S>
void mcmc_pack_edge_cut(context_handle_t<TenT>& ctx,TenT& tensor,
    const std::vector<std::pair<S,S>>& edges,const mcmc_edge_ports& ports,S site) {
    const auto cut=mcmc_edge_cut(edges,ports,site);const auto shape=tcapi::shape(ctx,tensor);
    if(shape.size()!=cut.size())throw std::invalid_argument("native cut rank mismatch");
    tcapi::List<tcapi::bond_idx_t<TenT>> permutation;tcapi::shape_t<TenT> packed;
    for(auto e:GetSurroundingBondIndex(site,edges)) {
        tcapi::bond_dim_t<TenT> dim=1;
        for(auto factor:ports[e]) {
            auto axis=std::lower_bound(cut.begin(),cut.end(),factor)-cut.begin();
            permutation.push_back(axis);
            if(shape[axis]<1 || dim>std::numeric_limits<tcapi::bond_dim_t<TenT>>::max()/shape[axis])
                throw std::overflow_error("native cut dimension");
            dim*=shape[axis];
        }
        packed.push_back(dim);
    }
    bool identity=true;for(std::size_t i=0;i<permutation.size();++i)identity=identity && permutation[i]==i;
    if(!identity)tcapi::transpose(ctx,tensor,permutation);
    if(packed!=shape)tcapi::reshape(ctx,tensor,packed);
}
template<class TenT,class S>
void mcmc_unpack_edge_cut(context_handle_t<TenT>& ctx,TenT& tensor,
    const std::vector<std::pair<S,S>>& edges,const mcmc_edge_ports& ports,S site,
    const std::map<std::size_t,unsigned long long>& dimensions) {
    std::vector<std::size_t> factors;tcapi::shape_t<TenT> expanded;
    for(auto e:GetSurroundingBondIndex(site,edges))for(auto f:ports[e]) {
        factors.push_back(f);expanded.push_back(dimensions.at(f));
    }
    if(tcapi::shape(ctx,tensor)!=expanded)tcapi::reshape(ctx,tensor,expanded);
    const auto sorted=mcmc_edge_cut(edges,ports,site);
    if(sorted!=factors) {
        tcapi::List<tcapi::bond_idx_t<TenT>> permutation;
        for(auto f:sorted)permutation.push_back(std::find(factors.begin(),factors.end(),f)-factors.begin());
        tcapi::transpose(ctx,tensor,permutation);
    }
}

/** Execute a native region, retaining common subtrees once and the target spine
 * once per state. Provider returns immutable canonical cuts on their owners.
 * Only operands and common neighbors are materialized, not the whole graph.
 * All steps and exterior-axis changes are compiled before numerical execution.
 */
template<class TenT,class S,class Provider>
std::vector<TenT> mcmc_execute_edge_region(context_handle_t<TenT>& ctx,
    const mcmc_edge_region_plan<S>& region,const std::map<S,S>& parent_owners,
    const std::vector<int>& participants,MPI_Comm comm,int count,Provider&& provider,
    int phase,std::vector<std::size_t>* audit=nullptr,
    const std::map<S,TenT*>* consume=nullptr) {
    if(count<1)throw std::invalid_argument("native candidate count");
    int rank;MPI_Comm_rank(comm,&rank);
    const auto& graph=std::get<0>(region);const auto& ports=std::get<1>(region);
    const auto& ops=std::get<2>(region);const auto target=std::get<4>(region);
    const auto& metadata=std::get<5>(region);const auto common=std::get<6>(region);
    auto owners=mcmc_region_local_owners(parent_owners,std::get<7>(region),participants);
    using Dimensions=std::map<std::size_t,unsigned long long>;
    std::map<S,Dimensions> dimensions;
    std::vector<TenT> base;std::vector<S> ids;
    for(auto site:std::get<7>(region))if(owners.at(site)==rank) {
        const auto& value=provider(site,0);const auto sh=tcapi::shape(ctx,value);
        const auto cut=mcmc_edge_cut(graph,ports,site);
        if(sh.size()!=cut.size())throw std::invalid_argument("native source cut mismatch");
        dimensions[site];
        for(std::size_t k=0;k<cut.size();++k)dimensions[site][cut[k]]=sh[k];
        // Caller grants a move only after considering later candidates and
        // windows. Raw caches and roulette-pending A_s never enter this map.
        if(consume && consume->count(site)) {
            base.push_back(std::move(*consume->at(site)));if(audit)++audit->at(11);
        }else {base.push_back(tcapi::copy(ctx,value));if(audit)++audit->at(7);}
        ids.push_back(site);mcmc_pack_edge_cut(ctx,base.back(),graph,ports,site);
    }
    // Factor shapes are communicated only between the same two operand owners,
    // once per step. Repeated candidates reuse this metadata. No line reduction.
    std::vector<Dimensions> result_dimensions(ops.size());
    auto execute=[&](std::size_t k,auto& work,auto& work_ids,auto& work_owners,auto& dims,bool first) {
        const auto st=ops[k];const auto& before=std::get<0>(metadata[k]);
        const auto& before_ports=std::get<1>(metadata[k]);const auto& plan=std::get<2>(metadata[k]);
        const auto& after_ports=std::get<3>(metadata[k]);
        const int ra=work_owners.at(st[0]),rb=work_owners.at(st[1]);
        const auto output=mcmc_edge_cut(plan.edges,after_ports,st[2]);
        if(first && !output.empty()) {
            const auto bcut=mcmc_edge_cut(before,before_ports,st[1]);
            std::vector<unsigned long long> bshape(bcut.size());
            if(rank==rb)for(std::size_t i=0;i<bcut.size();++i)bshape[i]=dims.at(st[1]).at(bcut[i]);
            if(ra!=rb) {
                if(rank==rb)MPI_Send(bshape.data(),static_cast<int>(bshape.size()),MPI_UNSIGNED_LONG_LONG,ra,107,comm);
                if(rank==ra)MPI_Recv(bshape.data(),static_cast<int>(bshape.size()),MPI_UNSIGNED_LONG_LONG,rb,107,comm,MPI_STATUS_IGNORE);
            }
            if(rank==ra) {
                auto merged=dims.at(st[0]);for(std::size_t i=0;i<bcut.size();++i)merged[bcut[i]]=bshape[i];
                for(auto f:output)result_dimensions[k][f]=merged.at(f);
            }
        }
        auto edges=before;
        graph_tensor_contraction(ctx,edges,work,work_ids,work_owners,st[0],st[1],st[2],comm);
        dims.erase(st[1]);if(rank==ra)dims[st[2]]=result_dimensions[k];
        if(audit) {++audit->at(phase);++audit->at(9);audit->at(3)=std::max(audit->at(3),work.size());}
    };
    for(std::size_t k=0;k<common;++k) {
        execute(k,base,ids,owners,dimensions,true);
        if(audit && phase==1)++audit->at(8);
    }
    const auto& repeated_graph=common<ops.size()?std::get<0>(metadata[common]):graph;
    const auto& repeated_ports=common<ops.size()?std::get<1>(metadata[common]):ports;
    std::vector<TenT> results;
    if(owners.at(target)==rank)results.resize(count);
    for(int state=0;state<count;++state) {
        std::vector<TenT> work;std::vector<S> work_ids;auto work_owners=owners;auto dims=dimensions;
        for(auto site:std::get<8>(region))if(owners.at(site)==rank) {
            auto it=std::find(ids.begin(),ids.end(),site);
            if(it==ids.end())throw std::logic_error("native retained input missing");
            const auto i=it-ids.begin();
            if(site==target && state>0) {
                work.push_back(tcapi::copy(ctx,provider(site,state)));
                mcmc_pack_edge_cut(ctx,work.back(),repeated_graph,repeated_ports,site);
                if(audit)++audit->at(7);
            } else if(site==target || state+1==count)work.push_back(std::move(base[i]));
            else {work.push_back(tcapi::copy(ctx,base[i]));if(audit)++audit->at(7);}
            work_ids.push_back(site);
        }
        for(std::size_t k=common;k<ops.size();++k)execute(k,work,work_ids,work_owners,dims,state==0);
        if(owners.at(target)==rank) {
            auto it=std::find(work_ids.begin(),work_ids.end(),target);
            auto& value=work.at(it-work_ids.begin());
            const auto& final_graph=ops.empty()?graph:std::get<2>(metadata.back()).edges;
            const auto& final_ports=ops.empty()?ports:std::get<3>(metadata.back());
            mcmc_unpack_edge_cut(ctx,value,final_graph,final_ports,target,dims.at(target));
            results[state]=std::move(value);
        }
    }
    if(audit)++audit->at(14); // Execution reused compiled metadata, including scalar joins.
    return results;
}

// Scalar joins have one composite edge, no exterior neighbors and no retained
// output cut. Prepare their local ownership and packed shapes once per window.
template<class TenT,class S>
using mcmc_edge_scalar_metadata=std::tuple<std::map<S,S>,std::vector<S>,
    std::vector<tcapi::shape_t<TenT>>,std::vector<int>>;

template<class TenT,class S>
mcmc_edge_scalar_metadata<TenT,S> mcmc_prepare_edge_scalar_join(
    context_handle_t<TenT>& ctx,const mcmc_edge_region_plan<S>& plan,
    const std::map<S,S>& parent_owners,const std::vector<int>& participants,
    int parent_rank,const TenT* a,const TenT* b) {
    const auto& edges=std::get<0>(plan);const auto& ports=std::get<1>(plan);
    const auto& steps=std::get<2>(plan);
    if(edges.size()!=1 || steps.size()!=1 || ports.size()!=1 ||
       !std::is_sorted(ports[0].begin(),ports[0].end()))
        throw std::invalid_argument("scalar join requires a canonical single-edge cut");
    const auto st=steps[0];
    auto owners=mcmc_region_local_owners(parent_owners,std::get<7>(plan),participants);
    std::vector<S> ids;std::vector<tcapi::shape_t<TenT>> packed;std::vector<int> operand;
    for(int side=0;side<2;++side)if(parent_owners.at(st[side])==parent_rank) {
        const auto* value=side==0?a:b;
        if(!value)throw std::invalid_argument("missing scalar join operand");
        const auto shape=tcapi::shape(ctx,*value);
        if(shape.size()!=ports[0].size())throw std::invalid_argument("scalar cut rank mismatch");
        tcapi::bond_dim_t<TenT> size=1;
        for(auto dim:shape) {
            if(dim<1 || size>std::numeric_limits<tcapi::bond_dim_t<TenT>>::max()/dim)
                throw std::overflow_error("scalar cut dimension");
            size*=dim;
        }
        ids.push_back(st[side]);packed.push_back({size});operand.push_back(side);
    }
    return {std::move(owners),std::move(ids),std::move(packed),std::move(operand)};
}

template<class TenT,class S>
std::vector<TenT> mcmc_execute_edge_scalar_join(context_handle_t<TenT>& ctx,
    const mcmc_edge_region_plan<S>& plan,const mcmc_edge_scalar_metadata<TenT,S>& metadata,
    const TenT* a,TenT* b,bool consume_b,MPI_Comm comm,std::vector<std::size_t>* audit=nullptr) {
    // The existing graph API consumes these working containers. Fixed maps and
    // axes are copied from the descriptor, not discovered again per candidate.
    auto edges=std::get<0>(plan);auto owners=std::get<0>(metadata);auto ids=std::get<1>(metadata);
    const auto& packed=std::get<2>(metadata);const auto& operand=std::get<3>(metadata);
    std::vector<TenT> work;work.reserve(ids.size());
    for(std::size_t i=0;i<ids.size();++i) {
        if(operand[i]==1 && consume_b) {
            work.push_back(std::move(*b));if(audit)++audit->at(11);
        } else {work.push_back(tcapi::copy(ctx,operand[i]==0?*a:*b));if(audit)++audit->at(7);}
        // All canonical factors of the single edge are already in order.
        // No permutation, factor-dimension map or output normalization is needed.
        if(std::get<1>(plan)[0].size()!=1)tcapi::reshape(ctx,work.back(),packed[i]);
    }
    const auto st=std::get<2>(plan)[0];
    graph_tensor_contraction(ctx,edges,work,ids,owners,st[0],st[1],st[2],comm);
    if(audit){++audit->at(1);++audit->at(9);++audit->at(14);audit->at(3)=std::max(audit->at(3),work.size());}
    return work;
}

/** Full-state heat bath on real quotient graphs of the original three-line
 * edges. Metadata is identical on all parent ranks; values stay on owners.
 */
template<class TenT,class S,class RNG>
void mcmc_update_line_heat_bath(context_handle_t<TenT>& ctx,
    const std::vector<std::pair<S,S>>& edges,const std::vector<TenT>& W,
    const std::vector<S>& site_idx,const std::map<S,S>& owners,
    const std::vector<TenT>& upper,const std::vector<S>& upper_ids,
    std::vector<TenT>& target,const std::vector<S>& target_ids,
    const std::vector<TenT>& lower,const std::vector<S>& lower_ids,
    const std::vector<S>& line,const std::vector<std::vector<S>>& groups,
    std::vector<S>& states,MPI_Comm parent,RNG& rng,int direction,int width,
    std::vector<std::size_t>* audit=nullptr,
    const std::map<std::vector<int>,MPI_Comm>* shared_communicators=nullptr,
    const mcmc_edge_line_plan<S>* prepared=nullptr) try {
    int rank;MPI_Comm_rank(parent,&rank);
    if(audit)audit->assign(15,0);
    mcmc_edge_line_plan<S> local_plan;
    if(!prepared) {local_plan=mcmc_build_edge_line_plan(edges,line,groups,owners,direction,width);prepared=&local_plan;}
    const auto& order=std::get<0>(*prepared);const auto& envplans=std::get<1>(*prepared);
    const auto& windows=std::get<2>(*prepared);
    if(windows.empty())return;
    std::map<std::vector<int>,MPI_Comm> temporary;
    if(!shared_communicators) {
        mcmc_acquire_region_communicators(parent,0,std::get<3>(*prepared),temporary);shared_communicators=&temporary;
    }
    bool participating=false;
    for(const auto& members:std::get<3>(*prepared))participating=participating || std::binary_search(members.begin(),members.end(),rank);
    if(!participating){mcmc_release_region_communicators(temporary);return;}
    std::vector<std::vector<TenT>> slices(line.size());std::vector<int> selected(line.size(),0);
    std::map<S,const TenT*> raw;
    for(std::size_t j=0;j<groups.size();++j)for(auto site:groups[j])if(owners.at(site)==rank) {
        const bool physical=site==line[j];const auto* source=&W;const auto* source_ids=&site_idx;
        if(!physical) {
            if(std::find(upper_ids.begin(),upper_ids.end(),site)!=upper_ids.end()){source=&upper;source_ids=&upper_ids;}
            else {source=&lower;source_ids=&lower_ids;}
        }
        const auto index=std::find(source_ids->begin(),source_ids->end(),site)-source_ids->begin();
        const auto& value=source->at(index);const auto sh=tcapi::shape(ctx,value);
        if(sh.size()!=GetSurroundingBondIndex(site,edges).size()+physical)throw std::invalid_argument("native raw tensor axes");
        if(!physical){raw[site]=&value;continue;}
        if(sh.back()<1 || sh.back()>std::numeric_limits<int>::max())throw std::invalid_argument("physical dimension");
        selected[j]=states.at(index);
        if(selected[j]<0 || selected[j]>=sh.back())throw std::invalid_argument("native state outside basis");
        if(static_cast<std::size_t>(std::find(target_ids.begin(),target_ids.end(),site)-target_ids.begin())>=target.size())throw std::invalid_argument("missing commit target");
        slices[j].resize(sh.back());
        tcapi::List<tcapi::bond_label_t<TenT>> labels(sh.size());std::iota(labels.begin(),labels.end(),0);
        for(std::size_t s=0;s<slices[j].size();++s) {
            tcapi::List<tcapi::bond_label_t<TenT>> out;
            mcmc_attach_window_state(ctx,value,labels,labels.back(),s,slices[j][s],out);
        }
    }
    auto raw_value=[&](S site)->const TenT& {
        auto it=std::find(line.begin(),line.end(),site);
        return it==line.end()?*raw.at(site):slices.at(it-line.begin()).at(0);
    };
    std::vector<std::vector<TenT>> suffix(line.size()+1);std::vector<TenT> prefix;
    if(audit){audit->at(2)=line.size();audit->at(5)=1;}
    for(std::size_t end=line.size();end>static_cast<std::size_t>(width);) {
        --end;const auto j=order[end];const auto& plan=envplans[end];const auto& participants=std::get<3>(plan);
        auto comm=shared_communicators->at(participants);if(comm==MPI_COMM_NULL)continue;
        if(audit){++audit->at(10);audit->at(12)+=std::get<2>(plan).size();}
        auto provider=[&](S site,int)->const TenT& {
            if(site==line[j])return slices[j].at(selected[j]);
            if(end+1<line.size() && site==line[order[end+1]])return suffix[end+1].at(0);
            return raw_value(site);
        };
        suffix[end]=mcmc_execute_edge_region<TenT>(ctx,plan,owners,participants,comm,1,provider,0,audit);
    }
    for(std::size_t pos=0;pos<windows.size();++pos) {
        const auto a=order[pos],b=order[pos+width-1],end=pos+width;
        const auto& win=windows[pos];const auto& participants=std::get<3>(win);
        auto comm=shared_communicators->at(participants);if(comm==MPI_COMM_NULL)continue;
        if(audit){audit->at(10)+=width;audit->at(12)+=std::get<2>(std::get<0>(win)).size()+std::get<2>(std::get<1>(win)).size();}
        auto local_rank=[&](S site){return static_cast<int>(std::lower_bound(participants.begin(),participants.end(),owners.at(site))-participants.begin());};
        int p=std::get<4>(win)[0],q=std::get<4>(win)[1];
        auto dimension=[&](int& d,std::size_t group) {
            if(d==0) {
                d=rank==owners.at(line[group])?static_cast<int>(slices[group].size()):0;
                MPI_Bcast(&d,1,MPI_INT,local_rank(line[group]),comm);
            } else if(rank==owners.at(line[group]) && d!=static_cast<int>(slices[group].size()))
                throw std::invalid_argument("physical dimension differs from prepared plan");
        };
        dimension(p,a);if(width==2)dimension(q,b);
        if(p<1 || q<1 || p>std::numeric_limits<int>::max()/q)throw std::overflow_error("window candidate count");
        auto aprovider=[&](S site,int state)->const TenT& {
            if(site==line[a])return slices[a].at(state);
            if(pos && site==line[order[pos-1]])return prefix.at(0);
            if(end<line.size() && site==line[order[end]])return suffix[end].at(0);
            return raw_value(site);
        };
        // This prefix version has no consumer after A's candidate preparation.
        std::map<S,TenT*> consume_a;
        if(pos && rank==owners.at(line[order[pos-1]]))consume_a[line[order[pos-1]]]=&prefix.at(0);
        auto candidates=mcmc_execute_edge_region<TenT>(ctx,std::get<0>(win),owners,participants,comm,p,aprovider,1,audit,&consume_a);
        std::vector<TenT> right;
        if(width==2) {
            auto bprovider=[&](S site,int state)->const TenT& {
                if(site==line[b])return slices[b].at(state);
                if(site==line[a])return candidates.at(0);
                if(end<line.size() && site==line[order[end]])return suffix[end].at(0);
                return raw_value(site);
            };
            // A's exterior-axis use has finished. This proper suffix will not
            // be read by a later window, so B may consume it now.
            std::map<S,TenT*> consume_b;
            if(end<line.size() && rank==owners.at(line[order[end]]))consume_b[line[order[end]]]=&suffix[end].at(0);
            right=mcmc_execute_edge_region<TenT>(ctx,std::get<1>(win),owners,participants,comm,q,bprovider,1,audit,&consume_b);
        }
        std::vector<tcapi::real_t<TenT>> amplitudes(rank==owners.at(line[a])?p*q:0),weights;int bad=0;
        const bool join=width==2 || end<line.size();
        mcmc_edge_scalar_metadata<TenT,S> scalar_metadata;
        if(join) {
            const auto partner=width==2?line[b]:line[order[end]];
            scalar_metadata=mcmc_prepare_edge_scalar_join(ctx,std::get<2>(win),owners,participants,rank,
                rank==owners.at(line[a])?&candidates.at(0):nullptr,
                rank==owners.at(partner)?(width==2?&right.at(0):&suffix[end].at(0)):nullptr);
            if(audit)++audit->at(13);
        }
        for(int t=0;t<q;++t)for(int s=0;s<p;++s) {
            std::vector<TenT> scalar;
            if(join) {
                const auto partner=width==2?line[b]:line[order[end]];
                scalar=mcmc_execute_edge_scalar_join(ctx,std::get<2>(win),scalar_metadata,
                    rank==owners.at(line[a])?&candidates.at(s):nullptr,
                    rank==owners.at(partner)?(width==2?&right.at(t):&suffix[end].at(0)):nullptr,
                    s+1==p,comm,audit);
            }
            if(rank==owners.at(line[a]))try {
                amplitudes[s+p*t]=std::abs(tcapi::get_elem(ctx,join?scalar.at(0):candidates.at(s),tcapi::elem_coors_t<TenT>{}));
                if(!std::isfinite(amplitudes[s+p*t]))bad=1;
            }catch(const std::exception&){bad=1;}
        }
        unsigned long long chosen=0;
        if(rank==owners.at(line[a]) && !bad)try {
            mcmc_heat_bath_weights(amplitudes,weights);chosen=mcmc_draw_heat_bath(weights,rng);
        }catch(const std::exception&){bad=1;}
        MPI_Bcast(&bad,1,MPI_INT,local_rank(line[a]),comm);
        if(bad)throw std::runtime_error("window scalar/weights failed");
        MPI_Bcast(&chosen,1,MPI_UNSIGNED_LONG_LONG,local_rank(line[a]),comm);
        selected[a]=chosen%p;if(width==2)selected[b]=chosen/p;
        prefix.clear();if(rank==owners.at(line[a]))prefix.push_back(std::move(candidates.at(selected[a])));
        if(audit)++audit->at(4);
        // Proper suffix[end] has no future window consumer; prefix is the only
        // selected leading candidate retained for the next overlapping window.
        suffix[end].clear();
    }
    mcmc_commit_window_state(slices,site_idx,owners,line,selected,rank,target,target_ids,states);
    mcmc_release_region_communicators(temporary);
} catch(const std::exception& e){detail::abort_tensor_operation(parent,e.what());}
  catch(...){detail::abort_tensor_operation(parent,"unknown native heat-bath failure");}
} // namespace tnbp
#endif
