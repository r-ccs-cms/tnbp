/// This file is a part of r-ccs-cms/tnbp
/** @file mcmc_heat_bath.h
 * @brief Distributed graph contractions for full-state heat-bath windows.
 */
#ifndef TNBP_PTNS_FUNCTION_MCMC_HEAT_BATH_H
#define TNBP_PTNS_FUNCTION_MCMC_HEAT_BATH_H
#include "tnbp/framework/graph_contraction.h"
#include "tnbp/ptns/function/mcmc_region_communicators.h"
#include "tnbp/framework/graph_region_contraction_plan.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
namespace tnbp {

/** @brief Fix one explicitly labeled physical axis, leaving sorted virtual axes. */
template<class TenT>
void mcmc_attach_window_state(context_handle_t<TenT>& ctx, const TenT& w,
    const tcapi::List<tcapi::bond_label_t<TenT>>& labels,
    tcapi::bond_label_t<TenT> physical, std::size_t state,
    TenT& out, tcapi::List<tcapi::bond_label_t<TenT>>& out_labels) {
    auto it=std::find(labels.begin(),labels.end(),physical);
    if(it==labels.end()) throw std::invalid_argument("missing heat-bath physical axis");
    auto d=tcapi::shape(ctx,w).at(std::distance(labels.begin(),it));
    if(state>=static_cast<std::size_t>(d)) throw std::invalid_argument("invalid heat-bath state");
    auto s=tcapi::zeros<TenT>(ctx,tcapi::shape_t<TenT>{d});
    tcapi::set_elem(ctx,s,tcapi::elem_coors_t<TenT>{static_cast<tcapi::bond_dim_t<TenT>>(state)},tcapi::elem_t<TenT>(1));
    out_labels=labels;out_labels.erase(out_labels.begin()+std::distance(labels.begin(),it));
    tcapi::contract(ctx,w,labels,s,tcapi::List<tcapi::bond_label_t<TenT>>{physical},out,out_labels);
}

/** @brief Stable full-state Born weights from finite absolute amplitudes.
 * A common scale preserves relative weights. All-zero amplitudes are an error;
 * no uniform fallback or charge-based threshold is applied.
 */
template<class RealT>
void mcmc_heat_bath_weights(const std::vector<RealT>& amplitudes, std::vector<RealT>& weights) {
    RealT scale=0;
    for(auto a:amplitudes) {
        if(!std::isfinite(a) || a<0) throw std::runtime_error("nonfinite/negative heat-bath amplitude");
        scale=std::max(scale,a);
    }
    if(scale==0) throw std::runtime_error("all heat-bath candidate amplitudes are zero");
    weights.resize(amplitudes.size());
    for(std::size_t k=0;k<amplitudes.size();++k) {
        auto r=amplitudes[k]/scale; weights[k]=r*r;
    }
}

/** @brief Roulette selection over all candidates, including the current state.
 * Uses extended-precision accumulation; the final positive bin covers rounding
 * at the endpoint. Zero bins are never selected. No tensor or MPI operations.
 */
template<class RealT, class RNG>
std::size_t mcmc_draw_heat_bath(const std::vector<RealT>& weights, RNG& rng) {
    RealT scale=0;
    for(auto w:weights) {
        if(!std::isfinite(w) || w<0) throw std::runtime_error("nonfinite/negative heat-bath weight");
        scale=std::max(scale,w);
    }
    if(scale==0) throw std::runtime_error("all heat-bath weights are zero");
    long double sum=0;
    std::size_t last=0;
    for(std::size_t k=0;k<weights.size();++k) {
        sum+=static_cast<long double>(weights[k])/scale;
        if(weights[k]>0) last=k;
    }
    const auto draw=std::generate_canonical<long double,std::numeric_limits<long double>::digits>(rng)*sum;
    long double cumulative=0;
    for(std::size_t k=0;k<weights.size();++k) {
        cumulative+=static_cast<long double>(weights[k])/scale;
        if(weights[k]>0 && draw<cumulative) return k;
    }
    return last;
}

/** @brief Move selected raw target slices and commit states on each owner.
 * Slices retain the original virtual-axis order and contain no environment.
 * No physical projection, candidate evaluation or contraction occurs here.
 */
template<class TenT, class IntT>
void mcmc_commit_window_state(std::vector<std::vector<TenT>>& slices,
    const std::vector<IntT>& site_idx,
    const std::map<IntT,IntT>& owners, const std::vector<IntT>& line,
    const std::vector<int>& selected, int rank,
    std::vector<TenT>& target, const std::vector<IntT>& target_ids,
    std::vector<IntT>& states) {
    for(std::size_t j=0;j<line.size();++j) if(owners.at(line[j])==rank) {
        auto wi=std::find(site_idx.begin(),site_idx.end(),line[j])-site_idx.begin();
        auto ti=std::find(target_ids.begin(),target_ids.end(),line[j])-target_ids.begin();
        target.at(ti)=std::move(slices.at(j).at(selected[j]));
        states.at(wi)=static_cast<IntT>(selected[j]);
    }
}


/** @brief Prepare a distributed three-line graph with explicit physical leaves.
 * W stays on its original owner. Each target physical axis is represented by
 * an appended graph edge to a metadata-only leaf; no implicit open axes
 * are passed to graph contraction. Boundary tensors retain their owners.
 * Graph, leaf IDs and dimensions are shared; numerical tensors are rank-local.
 */
template<class TenT,class IntT>
void mcmc_prepare_heat_bath_graph(context_handle_t<TenT>& ctx,
    const std::vector<std::pair<IntT,IntT>>& edges,const std::vector<TenT>& W,
    const std::vector<IntT>& site_idx,const std::map<IntT,IntT>& owners,
    const std::vector<TenT>& upper,const std::vector<IntT>& upper_ids,
    const std::vector<TenT>& lower,const std::vector<IntT>& lower_ids,
    const std::vector<IntT>& line,const std::vector<std::vector<IntT>>& groups,
    const std::vector<IntT>& states,MPI_Comm comm,
    std::vector<std::pair<IntT,IntT>>& graph,std::vector<TenT>& tensors,
    std::vector<IntT>& ids,std::map<IntT,IntT>& ranks,
    std::vector<IntT>& leaves,std::vector<int>& dimensions) try {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    std::set<IntT> members;
    if(groups.size()!=line.size()) throw std::invalid_argument("heat-bath group count");
    for(std::size_t j=0;j<groups.size();++j) {
        if(groups[j].empty() || groups[j].front()!=line[j]) throw std::invalid_argument("heat-bath group anchor");
        for(auto site:groups[j]) if(!members.insert(site).second)
            throw std::invalid_argument("duplicate heat-bath group member");
    }
    int bad=0;
    graph=edges;tensors.clear();ids.clear();ranks.clear();
    dimensions.assign(line.size(),0);leaves.clear();
    for(auto site:members) {
        auto owner=owners.find(site);
        if(owner==owners.end() || owner->second<0 || owner->second>=size) {bad=1;continue;}
        ranks[site]=owner->second;
        if(owner->second!=rank) continue;
        const auto target=std::find(line.begin(),line.end(),site);
        const bool is_target=target!=line.end();
        const auto* source=&W;const auto* source_ids=&site_idx;
        if(!is_target) {
            if(std::find(upper_ids.begin(),upper_ids.end(),site)!=upper_ids.end()) {source=&upper;source_ids=&upper_ids;}
            else {source=&lower;source_ids=&lower_ids;}
        }
        auto found=std::find(source_ids->begin(),source_ids->end(),site);
        auto index=static_cast<std::size_t>(found-source_ids->begin());
        if(found==source_ids->end() || index>=source->size()) {bad=1;continue;}
        auto shape=tcapi::shape(ctx,source->at(index));
        auto degree=GetSurroundingBondIndex(site,edges).size();
        if(shape.size()!=degree+(is_target ? 1:0)) {bad=1;continue;}
        if(is_target) {
            if(shape.back()<1 || shape.back()>std::numeric_limits<int>::max() || index>=states.size() ||
               states[index]<0 || states[index]>=shape.back()) {bad=1;continue;}
            dimensions[target-line.begin()]=static_cast<int>(shape.back());
        }
        tensors.push_back(tcapi::copy(ctx,source->at(index)));ids.push_back(site);
    }
    for(auto e:edges) if(!members.count(e.first) || !members.count(e.second)) bad=1;
    if(bad) throw std::invalid_argument("heat-bath graph/tensor/state/owner mismatch");
    MPI_Allreduce(MPI_IN_PLACE,dimensions.data(),static_cast<int>(dimensions.size()),MPI_INT,MPI_MAX,comm);
    std::set<IntT> occupied;for(auto entry:owners) occupied.insert(entry.first);
    IntT next=std::numeric_limits<IntT>::max();
    for(std::size_t j=0;j<line.size();++j) {
        while(occupied.count(next)) {
            if(next==std::numeric_limits<IntT>::lowest()) throw std::overflow_error("no physical leaf label available");
            --next;
        }
        const auto leaf=next;occupied.insert(leaf);leaves.push_back(leaf);
        graph.emplace_back(line[j],leaf);
    }
} catch(const std::exception& e) {detail::abort_tensor_operation(comm,e.what());}
  catch(...) {detail::abort_tensor_operation(comm,"unknown heat-bath failure");}


/** @brief Put cache axes in original edge-token order, without fusing cut legs. */
template<class TenT,class IntT>
void mcmc_canonicalize_cache(context_handle_t<TenT>& ctx,std::vector<TenT>& local,
                            std::vector<IntT>& axes) {
    auto sorted=axes;std::sort(sorted.begin(),sorted.end());
    if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())
        throw std::invalid_argument("duplicate cached boundary leg");
    if(sorted!=axes && !local.empty()) {
        tcapi::List<tcapi::bond_idx_t<TenT>> permutation;
        for(auto token:sorted) permutation.push_back(std::find(axes.begin(),axes.end(),token)-axes.begin());
        tcapi::transpose(ctx,local.at(0),permutation);
    }
    axes=std::move(sorted);
}

// Local failure is fatal; successful preparation/publication has no collective.
inline void mcmc_plan_check(MPI_Comm comm,int bad,const char* message) {
    if(bad) detail::abort_tensor_operation(comm,message);
}

/** Execute a prepared region join; only private operand copies are mutated.
 * External graph leaves have metadata but no numerical tensors or owners.
 */
template<class TenT,class TokenT>
void mcmc_contract_cached_tensors_prepared(context_handle_t<TenT>& ctx,
    const TenT* a,int owner_a,const TenT* b,int owner_b,
    const std::vector<std::pair<int,int>>& join_edges,
    const std::map<int,std::vector<tcapi::bond_dim_t<TenT>>>& join_shapes,
    const std::vector<std::pair<std::size_t,std::size_t>>& fuse_a,
    const std::vector<std::pair<std::size_t,std::size_t>>& fuse_b,
    const std::vector<TokenT>& result_axes,MPI_Comm comm,std::vector<TenT>& result,
    int phase=1,std::vector<std::size_t>* audit=nullptr,
    std::vector<std::vector<std::size_t>>* trace=nullptr) {
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);int bad=0;
    std::vector<std::pair<int,int>> graph;std::vector<TenT> tensors,out;
    std::vector<int> ids;std::map<int,int> owners;
    std::uint64_t macs=0,elements=0;
    try {
        if(owner_a<0 || owner_a>=size || owner_b<0 || owner_b>=size ||
           (a!=nullptr)!=(rank==owner_a) || (b!=nullptr)!=(rank==owner_b))
            throw std::invalid_argument("cached region is not owner-local");
        graph=join_edges;owners={{0,owner_a},{1,owner_b}};tensors.reserve(2);ids.reserve(2);out.reserve(1);
        auto append=[&](const TenT* input,int id,const auto& fusions) {
            if(!input) return;
            auto t=tcapi::copy(ctx,*input);
            for(auto axes:fusions) detail::fuse_graph_bond_axes(ctx,t,axes.first,axes.second);
            auto shape=tcapi::shape(ctx,t);
            if(std::vector<tcapi::bond_dim_t<TenT>>(shape.begin(),shape.end())!=join_shapes.at(id))
                throw std::invalid_argument("prepared region shape mismatch");
            tensors.push_back(std::move(t));ids.push_back(id);
        };
        append(a,0,fuse_a);append(b,1,fuse_b);
        if(trace) {
            auto e=join_edges;auto sh=join_shapes;
            macs=graph_tensor_contraction_cost(e,sh,0,1,0);elements=detail::graph_cost_elements(sh.at(0));
        }
    } catch(const std::exception&) {bad=1;}
    mcmc_plan_check(comm,bad,"cached join preparation failed");
    graph_tensor_contraction(ctx,graph,tensors,ids,owners,0,1,0,comm);
    try {
        if(rank==owner_a) {
            const auto index=std::find(ids.begin(),ids.end(),0)-ids.begin();
            if(static_cast<std::size_t>(index)>=tensors.size()) throw std::runtime_error("missing join result");
            const auto shape=tcapi::shape(ctx,tensors[index]);
            if(shape.size()!=result_axes.size()) throw std::runtime_error("join result rank");
            std::size_t axis=0;
            for(const auto& edge:join_edges) {
                if(edge.first==0 && edge.second==1) continue;
                if(shape.at(axis++)!=join_shapes.at(edge.second).at(0)) throw std::runtime_error("join result dimension");
            }
            if(trace && detail::graph_cost_elements(shape)!=elements)
                throw std::runtime_error("join element estimate mismatch");
            out.push_back(std::move(tensors[index]));
        }
        if(audit) {audit->at(phase)++;audit->at(3)=std::max(audit->at(3),std::size_t(2));audit->at(7)+=2;audit->at(9)++;}
        if(trace) trace->push_back({static_cast<std::size_t>(owner_a),static_cast<std::size_t>(owner_b),
            fuse_a.size()+1,static_cast<std::size_t>(macs),static_cast<std::size_t>(elements)});
    } catch(const std::exception&) {bad=1;}
    mcmc_plan_check(comm,bad,"cached join result preparation failed");
    result.swap(out);
}

/** Compatibility entry: prepare metadata once for this binary call. Group and
 * scalar candidate loops use the prepared entry directly to reuse metadata.
 */
template<class TenT,class IntT>
void mcmc_contract_cached_tensors(context_handle_t<TenT>& ctx,
    const TenT* a,const std::vector<IntT>& axes_a,int owner_a,
    const TenT* b,const std::vector<IntT>& axes_b,int owner_b,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<TenT>& result,std::vector<IntT>& axes_result,
    int phase=1,std::vector<std::size_t>* audit=nullptr,bool invalid_input=false,
    std::vector<std::vector<std::size_t>>* trace=nullptr) {
    std::vector<std::pair<int,int>> graph;std::map<int,std::vector<tcapi::bond_dim_t<TenT>>> shapes;
    std::vector<std::pair<std::size_t,std::size_t>> fa,fb;std::vector<IntT> axes;int bad=invalid_input;
    try {make_region_join_metadata(axes_a,axes_b,dimensions,graph,shapes,fa,fb,axes);}
    catch(const std::exception&) {bad=1;}
    mcmc_plan_check(comm,bad,"cached join metadata failed");
    std::vector<TenT> out;
    mcmc_contract_cached_tensors_prepared(ctx,a,owner_a,b,owner_b,graph,shapes,fa,fb,axes,comm,out,phase,audit,trace);
    result.swap(out);axes_result.swap(axes);
}

/** @brief Borrow owner-local region tensors; clone each only inside the graph join. */
template<class TenT,class IntT>
void mcmc_contract_cached_regions(context_handle_t<TenT>& ctx,
    const std::vector<TenT>& a,const std::vector<IntT>& axes_a,int owner_a,
    const std::vector<TenT>& b,const std::vector<IntT>& axes_b,int owner_b,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<TenT>& result,std::vector<IntT>& axes_result,
    int phase=1,std::vector<std::size_t>* audit=nullptr,
    std::vector<std::vector<std::size_t>>* trace=nullptr) {
    int rank;MPI_Comm_rank(comm,&rank);
    int invalid=a.size()!=static_cast<std::size_t>(rank==owner_a) || b.size()!=static_cast<std::size_t>(rank==owner_b);
    mcmc_contract_cached_tensors(ctx,a.empty() ? nullptr:&a[0],axes_a,owner_a,
        b.empty() ? nullptr:&b[0],axes_b,owner_b,dimensions,comm,result,axes_result,phase,audit,invalid!=0,trace);
}

/** @brief Slice only raw targets; retain boundary members separately on owners.
 * No complete group tensor is formed. Original edge tokens label every axis.
 */
template<class TenT,class IntT>
void mcmc_prepare_heat_bath_cache(context_handle_t<TenT>& ctx,
    const std::vector<std::pair<IntT,IntT>>& edges,const std::vector<TenT>& W,
    const std::vector<IntT>& site_idx,const std::map<IntT,IntT>& owners,
    const std::vector<TenT>& upper,const std::vector<IntT>& upper_ids,
    const std::vector<TenT>& lower,const std::vector<IntT>& lower_ids,
    const std::vector<IntT>& line,const std::vector<std::vector<IntT>>& groups,
    const std::vector<IntT>& states,MPI_Comm comm,
    std::vector<std::vector<TenT>>& slices,std::vector<std::vector<IntT>>& axes,
    std::vector<int>& physical_dimensions,std::vector<tcapi::bond_dim_t<TenT>>& bond_dimensions,
    std::vector<std::size_t>* audit=nullptr,
    std::vector<std::vector<std::vector<TenT>>>* boundary=nullptr,
    std::vector<std::vector<std::vector<IntT>>>* boundary_axes=nullptr,
    std::vector<std::vector<int>>* boundary_owners=nullptr) try {
    int rank;MPI_Comm_rank(comm,&rank);
    std::vector<std::pair<IntT,IntT>> graph;std::vector<TenT> raw;
    std::vector<IntT> ids,leaves;std::map<IntT,IntT> ranks;
    mcmc_prepare_heat_bath_graph(ctx,edges,W,site_idx,owners,upper,upper_ids,lower,lower_ids,
        line,groups,states,comm,graph,raw,ids,ranks,leaves,physical_dimensions);
    std::vector<IntT> required;
    for(const auto& group:groups) required.insert(required.end(),group.begin(),group.end());
    detail::validate_graph_tensors(ctx,required,graph,raw,ids,ranks,comm,false,false);
    std::vector<unsigned long long> dimensions(graph.size(),0);
    for(std::size_t j=0;j<ids.size();++j) {
        auto bonds=GetSurroundingBondIndex(ids[j],graph);auto shape=tcapi::shape(ctx,raw[j]);
        for(std::size_t k=0;k<bonds.size();++k) dimensions[bonds[k]]=static_cast<unsigned long long>(shape[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE,dimensions.data(),static_cast<int>(dimensions.size()),MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
    // Reuse the shape sharing required by the cache. Every endpoint must equal
    // the shared maximum; a smaller local dimension detects disagreement without
    // another reduction, including when the two endpoints have different owners.
    for(std::size_t j=0;j<ids.size();++j) {
        const auto bonds=GetSurroundingBondIndex(ids[j],graph);
        const auto shape=tcapi::shape(ctx,raw[j]);
        for(std::size_t k=0;k<bonds.size();++k)
            if(static_cast<unsigned long long>(shape[k])!=dimensions[bonds[k]])
                throw std::invalid_argument("heat-bath bond dimensions disagree with shared shape");
    }
    bond_dimensions.assign(dimensions.begin(),dimensions.end());
    if(boundary) {boundary->clear();boundary->resize(line.size());}
    if(boundary_axes) boundary_axes->assign(line.size(),{});
    if(boundary_owners) boundary_owners->assign(line.size(),{});
    slices.clear();slices.resize(line.size());axes.clear();axes.resize(line.size());
    for(std::size_t j=0;j<line.size();++j) {
        if(groups[j].size()>1 && (!boundary || !boundary_axes || !boundary_owners))
            throw std::invalid_argument("boundary cache outputs required");
        if(audit) audit->at(2)++;
        for(std::size_t k=0;k<groups[j].size();++k) {
            const auto site=groups[j][k];
            std::vector<IntT> tokens;
            for(auto e:GetSurroundingBondIndex(site,edges)) tokens.push_back(static_cast<IntT>(e));
            if(k==0) {
                axes[j]=tokens;
                if(rank==owners.at(site)) {
                    const auto index=std::find(ids.begin(),ids.end(),site)-ids.begin();
                    const auto physical=static_cast<IntT>(find_edge(graph,site,leaves[j])-graph.begin());
                    auto labels=tokens;labels.push_back(physical);
                    tcapi::List<tcapi::bond_label_t<TenT>> in(labels.begin(),labels.end()),out;
                    slices[j].resize(physical_dimensions[j]);
                    for(int state=0;state<physical_dimensions[j];++state)
                        mcmc_attach_window_state(ctx,raw.at(index),in,physical,state,slices[j][state],out);
                }
            } else {
                boundary->at(j).emplace_back();boundary_axes->at(j).push_back(std::move(tokens));
                boundary_owners->at(j).push_back(static_cast<int>(owners.at(site)));
                if(rank==owners.at(site)) {
                    const auto index=std::find(ids.begin(),ids.end(),site)-ids.begin();
                    boundary->at(j).back().push_back(std::move(raw.at(index)));
                }
            }
        }
    }
} catch(const std::exception& e) {detail::abort_tensor_operation(comm,e.what());}
  catch(...) {detail::abort_tensor_operation(comm,"unknown heat-bath failure");}

/** Build the old connected absorption tree without touching tensor values. */
template<class TokenT>
std::vector<std::array<std::size_t,3>> mcmc_group_fallback(
    const std::vector<std::vector<TokenT>>& leaves,bool has_side) {
    auto axes=leaves;const auto n=leaves.size();std::vector<std::array<std::size_t,3>> ops;
    const auto begin=std::size_t(1)+has_side;std::vector<bool> used(n,false);
    bool have=has_side;std::size_t current=has_side?1:0;
    auto join=[&](std::size_t a,std::size_t b) {
        std::vector<TokenT> out;
        std::set_symmetric_difference(axes[a].begin(),axes[a].end(),axes[b].begin(),axes[b].end(),std::back_inserter(out));
        const auto c=axes.size();ops.push_back({a,b,c});axes.push_back(std::move(out));return c;
    };
    auto attach=[&](std::size_t k) {current=have?join(k,current):k;have=true;used[k]=true;};
    for(std::size_t k=begin;k<n;++k)
        if(detail::regions_connected(leaves[has_side?1:0],leaves[k])) {attach(k);break;}
    while(have && !detail::regions_connected(axes[current],axes[0])) {
        std::size_t k=begin;for(;k<n;++k) if(!used[k] && detail::regions_connected(axes[current],axes[k])) break;
        if(k==n) throw std::invalid_argument("group cannot reach target");attach(k);
    }
    current=have?join(0,current):0;
    for(std::size_t remaining=std::count(used.begin()+begin,used.end(),false);remaining;--remaining) {
        std::size_t k=begin;for(;k<n;++k) if(!used[k] && detail::regions_connected(axes[current],axes[k])) break;
        if(k==n) throw std::invalid_argument("disconnected group");used[k]=true;current=join(current,k);
    }
    return ops;
}

template<class TokenT,class DimT>
std::uint64_t mcmc_prepare_group_contraction_plan(
    const std::vector<TokenT>& target_axes,int target_owner,
    const std::vector<std::vector<TokenT>>& boundary_axes,const std::vector<int>& boundary_owners,
    const std::vector<TokenT>& side_axes,bool has_side,int side_owner,
    const std::vector<DimT>& dimensions,int count,
    std::vector<std::vector<TokenT>>& leaf_axes,std::vector<int>& leaf_owners,
    std::vector<std::array<std::size_t,3>>& steps,std::uint64_t& max_elements,bool& used_fallback) {
    if(count<1 || boundary_axes.size()!=boundary_owners.size()) throw std::invalid_argument("group plan inputs");
    std::vector<std::vector<TokenT>> a{target_axes};std::vector<int> o{target_owner};
    if(has_side) {a.push_back(side_axes);o.push_back(side_owner);}
    a.insert(a.end(),boundary_axes.begin(),boundary_axes.end());o.insert(o.end(),boundary_owners.begin(),boundary_owners.end());
    detail::validate_region_inputs(a,dimensions,o,0);
    auto baseline=mcmc_group_fallback(a,has_side);std::vector<std::array<std::size_t,3>> ops;
    std::uint64_t peak=0;bool fallback=false;
    const auto cost=make_small_region_contraction_plan(a,dimensions,o,0,count,baseline,ops,peak,fallback);
    leaf_axes.swap(a);leaf_owners.swap(o);steps.swap(ops);max_elements=peak;used_fallback=fallback;return cost;
}

/** Execute common subtrees once, then the target-dependent spine per state.
 * Input tensors are borrowed. Public outputs are swapped only after local
 * success; any preparation or contraction error aborts the job. No environment is reconstructed.
 */
template<class TenT,class TokenT>
void mcmc_execute_group_contraction_plan(context_handle_t<TenT>& ctx,
    const std::vector<const std::vector<TenT>*>& inputs,int count,int state_offset,
    const std::vector<std::array<std::size_t,3>>& ops,const std::vector<std::vector<TokenT>>& axes,
    const std::vector<int>& owners,const std::vector<unsigned char>& dependent,
    const std::vector<std::size_t>& parent,const std::vector<std::size_t>& frontier,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<std::vector<TenT>>& candidates,std::vector<TokenT>& candidate_axes,
    int phase=1,std::vector<std::size_t>* audit=nullptr,std::vector<int>* steps=nullptr,
    std::vector<std::vector<std::size_t>>* trace=nullptr) {
    using Dim=tcapi::bond_dim_t<TenT>;
    int rank,size;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);int bad=0;
    const auto n=inputs.size();
    std::vector<std::vector<std::pair<int,int>>> graphs;
    std::vector<std::map<int,std::vector<Dim>>> shapes;
    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> fa,fb;
    std::vector<std::vector<TenT>> storage,work;std::vector<TokenT> out_axes;
    std::vector<unsigned char> retained;
    try {
        if(!n || count<1 || state_offset<0 || ops.size()!=n-1 || axes.size()!=2*n-1 ||
           owners.size()!=axes.size() || dependent.size()!=axes.size() || parent.size()!=axes.size())
            throw std::invalid_argument("invalid group execution metadata");
        storage.resize(axes.size());work.resize(count);out_axes=axes.back();retained.assign(axes.size(),0);
        if(steps) {
            const auto extra=detail::graph_cost_product(std::max(std::size_t(1),ops.size()),static_cast<std::uint64_t>(count));
            if(extra>steps->max_size()-steps->size()) throw std::length_error("group trace too large");
            steps->reserve(steps->size()+static_cast<std::size_t>(extra));
        }
        for(auto id:frontier) retained.at(id)=1;
        for(std::size_t i=0;i<n;++i) {
            if(owners[i]<0 || owners[i]>=size || (inputs[i]!=nullptr)!=(rank==owners[i]))
                throw std::invalid_argument("group input owner");
            if(inputs[i]) {
                const auto required=i==0?static_cast<std::size_t>(state_offset)+static_cast<std::size_t>(count):1;
                if(inputs[i]->size()<required || (i && inputs[i]->size()!=1)) throw std::invalid_argument("missing group input");
            }
        }
        graphs.resize(ops.size());shapes.resize(ops.size());fa.resize(ops.size());fb.resize(ops.size());
        for(std::size_t k=0;k<ops.size();++k) {
            const auto op=ops[k];std::vector<TokenT> result_axes;
            make_region_join_metadata(axes.at(op[0]),axes.at(op[1]),dimensions,graphs[k],shapes[k],fa[k],fb[k],result_axes);
            if(result_axes!=axes.at(op[2]) || owners.at(op[2])!=owners.at(op[0])) throw std::invalid_argument("group plan mismatch");
        }
        if(audit && audit->size()>12) audit->at(12)+=ops.size();
    } catch(const std::exception&) {bad=1;}
    mcmc_plan_check(comm,bad,"group execution preparation failed");
    auto input=[&](std::size_t id,int state)->const TenT* {
        if(rank!=owners[id]) return nullptr;
        if(id<n) return &inputs[id]->at(id==0?static_cast<std::size_t>(state_offset)+state:0);
        return &storage[id].at(0);
    };
    auto execute=[&](std::size_t k,int state) {
        const auto op=ops[k];const TenT *a=nullptr,*b=nullptr;
        // Shape/owner availability follows the validated tree and previous joins.
        a=input(op[0],state);b=input(op[1],state);
        mcmc_contract_cached_tensors_prepared(ctx,a,owners[op[0]],b,owners[op[1]],
            graphs[k],shapes[k],fa[k],fb[k],axes[op[2]],comm,storage[op[2]],phase,audit,trace);
        // No allocation here: diagnostic capacity was reserved in preflight.
        for(auto id:{op[0],op[1]}) if(id>=n && !retained[id]) storage[id].clear();
        if(steps) steps->push_back(dependent[op[2]]?-1:0);
    };
    for(std::size_t k=0;k<ops.size();++k) if(!dependent[ops[k][2]]) execute(k,0);
    for(int state=0;state<count;++state) {
        if(n==1) {
            try {
                if(rank==owners[0]) {
                    auto t=tcapi::copy(ctx,*input(0,state));
                    auto sh=tcapi::shape(ctx,t);
                    if(sh.size()!=axes[0].size()) throw std::invalid_argument("singleton shape rank");
                    for(std::size_t j=0;j<sh.size();++j) if(sh[j]!=dimensions.at(axes[0][j])) throw std::invalid_argument("singleton shape dimension");
                    work[state].push_back(std::move(t));
                }
                if(audit && phase==1) audit->at(8)++;
                if(steps) steps->push_back(-1);
            } catch(const std::exception&) {bad=1;break;}
        } else {
            for(std::size_t k=0;k<ops.size();++k) if(dependent[ops[k][2]]) execute(k,state);
            work[state].swap(storage.back());
        }
    }
    // Publish only completed local results, including the singleton path.
    mcmc_plan_check(comm,bad,"group candidate publication failed");
    candidates.swap(work);candidate_axes.swap(out_axes);
}

// Plain arrays stored together; no tensor data is retained in the plan cache.
template<class TokenT>
using mcmc_region_plan_arrays = std::tuple<
    std::vector<std::array<std::size_t,3>>, std::vector<std::vector<TokenT>>,
    std::vector<int>, std::vector<unsigned char>, std::vector<std::size_t>,
    std::vector<std::size_t>, std::uint64_t>;

/** Plan a group once and reuse its common subtrees across all target states.
 * Optional line-local key/plan arrays reuse searches for identical metadata.
 * steps is a diagnostic stage trace: 0 for common, -1 for target-dependent joins.
 */
template<class TenT,class IntT>
void mcmc_prepare_group_candidates(context_handle_t<TenT>& ctx,
    const std::vector<TenT>& target,const std::vector<IntT>& target_axes,int target_owner,int count,
    const std::vector<std::vector<TenT>>& boundary,const std::vector<std::vector<IntT>>& boundary_axes,
    const std::vector<int>& boundary_owners,
    const std::vector<TenT>& side,const std::vector<IntT>& side_axes,bool has_side,int side_owner,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<std::vector<TenT>>& candidates,std::vector<IntT>& candidate_axes,
    int phase=1,std::vector<std::size_t>* audit=nullptr,std::vector<int>* steps=nullptr,
    int state_offset=0,std::vector<std::vector<std::size_t>>* trace=nullptr,
    std::vector<std::vector<std::uint64_t>>* plan_keys=nullptr,
    std::vector<mcmc_region_plan_arrays<IntT>>* plan_cache=nullptr,
    std::uint64_t* total_macs=nullptr) {
    int rank;MPI_Comm_rank(comm,&rank);int bad=0;
    const mcmc_region_plan_arrays<IntT>* cached=nullptr;std::uint64_t total=0;
    std::vector<std::vector<IntT>> leaves,axes;std::vector<int> leaf_owners,owners;
    std::vector<unsigned char> dependent;std::vector<std::size_t> parent,frontier;
    std::vector<std::uint64_t> macs,elements,key;
    std::vector<std::array<std::size_t,3>> ops;std::vector<const std::vector<TenT>*> inputs;
    try {
        if(count<1 || boundary.size()!=boundary_axes.size() || boundary.size()!=boundary_owners.size() ||
           (plan_keys==nullptr)!=(plan_cache==nullptr)) throw std::invalid_argument("group inputs");
        leaves.push_back(target_axes);leaf_owners.push_back(target_owner);inputs.push_back(rank==target_owner?&target:nullptr);
        if(has_side) {leaves.push_back(side_axes);leaf_owners.push_back(side_owner);inputs.push_back(rank==side_owner?&side:nullptr);}
        for(std::size_t i=0;i<boundary.size();++i) {
            leaves.push_back(boundary_axes[i]);leaf_owners.push_back(boundary_owners[i]);
            inputs.push_back(rank==boundary_owners[i]?&boundary[i]:nullptr);
        }
        detail::validate_region_inputs(leaves,dimensions,leaf_owners,0);
        key={static_cast<std::uint64_t>(count),static_cast<std::uint64_t>(has_side),leaves.size()};
        for(std::size_t i=0;i<leaves.size();++i) {
            key.push_back(leaf_owners[i]);key.push_back(leaves[i].size());
            for(auto t:leaves[i]) {key.push_back(t);key.push_back(dimensions.at(t));}
        }
        std::size_t index=0;
        if(plan_keys) {
            if(plan_keys->size()!=plan_cache->size()) throw std::invalid_argument("plan cache mismatch");
            index=std::find(plan_keys->begin(),plan_keys->end(),key)-plan_keys->begin();
        }
        if(plan_keys && index<plan_keys->size()) {
            cached=&plan_cache->at(index);total=std::get<6>(*cached);if(audit && audit->size()>14) audit->at(14)++;
        } else {
            std::uint64_t peak=0;bool fallback=false;
            total=mcmc_prepare_group_contraction_plan(target_axes,target_owner,boundary_axes,boundary_owners,
                side_axes,has_side,side_owner,dimensions,count,leaves,leaf_owners,ops,peak,fallback);
            if(audit && audit->size()>11) {audit->at(10)++;audit->at(11)+=fallback;}
            describe_region_contraction_plan(leaves,dimensions,leaf_owners,0,ops,axes,owners,dependent,parent,frontier,macs,elements);
            if(plan_keys) {
                auto saved=std::make_tuple(std::move(ops),std::move(axes),std::move(owners),
                    std::move(dependent),std::move(parent),std::move(frontier),total);
                plan_keys->reserve(plan_keys->size()+1);plan_cache->reserve(plan_cache->size()+1);
                plan_keys->push_back(std::move(key));plan_cache->push_back(std::move(saved));cached=&plan_cache->back();
            }
        }
    } catch(const std::exception&) {bad=1;}
    mcmc_plan_check(comm,bad,"group contraction planning failed");
    mcmc_execute_group_contraction_plan(ctx,inputs,count,state_offset,
        cached?std::get<0>(*cached):ops,cached?std::get<1>(*cached):axes,
        cached?std::get<2>(*cached):owners,cached?std::get<3>(*cached):dependent,
        cached?std::get<4>(*cached):parent,cached?std::get<5>(*cached):frontier,
        dimensions,comm,candidates,candidate_axes,phase,audit,steps,trace);
    if(total_macs) *total_macs=total;
}

/** @brief Build all proper opposite-side environments once for this sweep.
 * order follows the sweep, so suffixes represent the right side for a forward
 * sweep and the left side for a reverse sweep. They contain only initial states
 * of untouched groups. No candidate evaluation calls this routine.
 */
template<class TenT,class IntT>
void mcmc_build_heat_bath_environments(context_handle_t<TenT>& ctx,
    const std::vector<std::vector<TenT>>& slices,const std::vector<std::vector<IntT>>& axes,
    const std::vector<int>& owners,const std::vector<int>& states,const std::vector<std::size_t>& order,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<std::vector<TenT>>& suffix,std::vector<std::vector<IntT>>& suffix_axes,
    int width=1,std::vector<std::size_t>* audit=nullptr,
    const std::vector<std::vector<std::vector<TenT>>>* boundary=nullptr,
    const std::vector<std::vector<std::vector<IntT>>>* boundary_axes=nullptr,
    const std::vector<std::vector<int>>* boundary_owners=nullptr,
    std::vector<std::vector<std::uint64_t>>* plan_keys=nullptr,
    std::vector<mcmc_region_plan_arrays<IntT>>* plan_cache=nullptr) {
    const auto n=order.size();
    if(width!=1 && width!=2) throw std::invalid_argument("environment window size");
    suffix.clear();suffix.resize(n+1);suffix_axes.clear();suffix_axes.resize(n+1);
    if(audit) audit->at(5)++;
    for(std::size_t pos=n;pos>static_cast<std::size_t>(width);--pos) {
        auto j=order[pos-1];
        const std::vector<std::vector<TenT>> empty_boundary;
        const std::vector<std::vector<IntT>> empty_axes;const std::vector<int> empty_owners;
        std::vector<std::vector<TenT>> candidate;
        mcmc_prepare_group_candidates(ctx,slices[j],axes[j],owners[j],1,
            boundary ? boundary->at(j):empty_boundary,boundary_axes ? boundary_axes->at(j):empty_axes,
            boundary_owners ? boundary_owners->at(j):empty_owners,
            suffix[pos],suffix_axes[pos],pos<n,pos<n ? owners[order[pos]]:owners[j],
            dimensions,comm,candidate,suffix_axes[pos-1],0,audit,nullptr,states[j],nullptr,plan_keys,plan_cache);
        suffix[pos-1]=std::move(candidate[0]);
    }
}

/** @brief Evaluate only the update window between two immutable cached sides.
 * No raw graph, raw group list, remaining line, or exterior contraction schedule
 * is accepted here. Cached region joins are distributed graph contractions.
 * Window boundary members remain separate until attached around each target.
 * Return all d_first*d_second weights plus reusable first-group intermediates.
 */
template<class TenT,class IntT>
void mcmc_evaluate_window_weights(context_handle_t<TenT>& ctx,
    const std::vector<TenT>& first,const std::vector<IntT>& first_axes,int first_owner,int first_dim,
    const std::vector<TenT>& second,const std::vector<IntT>& second_axes,int second_owner,int second_dim,bool pair,
    const std::vector<TenT>& prefix,const std::vector<IntT>& prefix_axes,bool has_prefix,int prefix_owner,
    const std::vector<TenT>& suffix,const std::vector<IntT>& suffix_axes,bool has_suffix,int suffix_owner,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<tcapi::real_t<TenT>>& amplitudes,std::vector<tcapi::real_t<TenT>>& weights,
    std::vector<std::vector<TenT>>& candidates,std::vector<IntT>& candidate_axes,
    std::vector<std::size_t>* audit=nullptr,
    const std::vector<std::vector<TenT>>& first_boundary={},
    const std::vector<std::vector<IntT>>& first_boundary_axes={},const std::vector<int>& first_boundary_owners={},
    const std::vector<std::vector<TenT>>& second_boundary={},
    const std::vector<std::vector<IntT>>& second_boundary_axes={},const std::vector<int>& second_boundary_owners={},
    std::vector<std::vector<std::uint64_t>>* plan_keys=nullptr,
    std::vector<mcmc_region_plan_arrays<IntT>>* plan_cache=nullptr) {
    using Real=tcapi::real_t<TenT>;int rank;MPI_Comm_rank(comm,&rank);
    if(first_dim<1 || second_dim<1 || static_cast<std::size_t>(first_dim)>
       static_cast<std::size_t>(std::numeric_limits<int>::max())/static_cast<std::size_t>(second_dim))
        throw std::overflow_error("heat-bath candidate count");
    std::uint64_t first_cost=0,second_cost=0;
    mcmc_prepare_group_candidates(ctx,first,first_axes,first_owner,first_dim,
        first_boundary,first_boundary_axes,first_boundary_owners,prefix,prefix_axes,has_prefix,prefix_owner,
        dimensions,comm,candidates,candidate_axes,1,audit,nullptr,0,nullptr,plan_keys,plan_cache,&first_cost);
    std::vector<std::vector<TenT>> right;std::vector<IntT> right_axes;
    if(pair) mcmc_prepare_group_candidates(ctx,second,second_axes,second_owner,second_dim,
        second_boundary,second_boundary_axes,second_boundary_owners,suffix,suffix_axes,has_suffix,suffix_owner,
        dimensions,comm,right,right_axes,1,audit,nullptr,0,nullptr,plan_keys,plan_cache,&second_cost);
    std::vector<std::pair<int,int>> graph;std::map<int,std::vector<tcapi::bond_dim_t<TenT>>> shapes;
    std::vector<std::pair<std::size_t,std::size_t>> fa,fb;std::vector<IntT> result_axes;
    std::vector<Real> work_amplitudes,work_weights;int bad=0;
    try {
        const auto group_cost=detail::region_cost_sum(first_cost,second_cost);
        work_amplitudes.assign(static_cast<std::size_t>(first_dim)*second_dim,Real(0));
        work_weights.resize(work_amplitudes.size());
        if(!pair && has_suffix) {
            int size;MPI_Comm_size(comm,&size);
            if(suffix_owner<0 || suffix_owner>=size || suffix.size()!=static_cast<std::size_t>(rank==suffix_owner))
                throw std::invalid_argument("window suffix is not owner-local");
        }
        if(pair || has_suffix) {
            make_region_join_metadata(candidate_axes,pair?right_axes:suffix_axes,dimensions,graph,shapes,fa,fb,result_axes);
            if(!result_axes.empty()) throw std::invalid_argument("window does not close");
            auto e=graph;auto sh=shapes;
            const auto cost=graph_tensor_contraction_cost(e,sh,0,1,0);
            detail::region_cost_sum(group_cost,detail::graph_cost_product(cost,static_cast<std::uint64_t>(first_dim)*second_dim));
            if(!sh.at(0).empty()) throw std::invalid_argument("window result is not scalar");
            if(audit && audit->size()>13) audit->at(13)++;
        } else if(!candidate_axes.empty()) throw std::invalid_argument("open singleton window");
    } catch(const std::exception&) {bad=1;}
    mcmc_plan_check(comm,bad,"window join preparation failed");
    for(int t=0;t<second_dim;++t) for(int state=0;state<first_dim;++state) {
        std::vector<TenT> evaluated;
        if(pair || has_suffix) {
            const auto& rhs=pair?right[t]:suffix;const int rhs_owner=pair?second_owner:suffix_owner;
            mcmc_contract_cached_tensors_prepared(ctx,rank==first_owner?&candidates[state].at(0):nullptr,first_owner,
                rank==rhs_owner?&rhs.at(0):nullptr,rhs_owner,graph,shapes,fa,fb,result_axes,comm,evaluated,1,audit);
        }
        if(rank==first_owner) try {
            const auto a=std::abs(tcapi::get_elem(ctx,(pair || has_suffix)?evaluated.at(0):candidates[state].at(0),tcapi::elem_coors_t<TenT>{}));
            if(!std::isfinite(a)) bad=1;
            work_amplitudes[static_cast<std::size_t>(state)+static_cast<std::size_t>(first_dim)*t]=a;
        } catch(const std::exception&) {bad=1;}
    }
    // Scalar reads can fail locally without affecting subsequent joins. Keep
    // completing the candidate loop, then notify once before sharing/drawing.
    if(rank==first_owner && !bad) try {mcmc_heat_bath_weights(work_amplitudes,work_weights);}
    catch(const std::exception&) {bad=1;}
    int failed=0;MPI_Bcast(&bad,1,MPI_INT,first_owner,comm);failed=bad;
    if(failed) throw std::runtime_error("window scalar/weights failed");
    MPI_Bcast(work_amplitudes.data(),static_cast<int>(work_amplitudes.size()),GetMpiType<Real>(),first_owner,comm);
    // Storage was reserved before the loop; owner already validated every value.
    if(rank!=first_owner) mcmc_heat_bath_weights(work_amplitudes,work_weights);
    amplitudes.swap(work_amplitudes);weights.swap(work_weights);
}

/** @brief Move the selected leading intermediate into the advancing environment.
 * Cached cut legs are unfused, so no distant neighbor arrays need copying.
 */
template<class TenT,class IntT>
void mcmc_absorb_group(std::size_t selected,std::vector<std::vector<TenT>>& candidates,
    const std::vector<IntT>& candidate_axes,std::vector<TenT>& prefix,std::vector<IntT>& prefix_axes,
    std::vector<std::size_t>* audit=nullptr) {
    prefix=std::move(candidates.at(selected));prefix_axes=candidate_axes;
    if(audit) audit->at(4)++;
}

} // namespace tnbp
#include "tnbp/ptns/function/mcmc_edge_heat_bath.h"
#endif
