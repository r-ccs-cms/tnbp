/// This file is a part of r-ccs-cms/tnbp
/** @file mcmc_heat_bath.h
 * @brief Distributed graph contractions for full-state heat-bath windows.
 */
#ifndef TNBP_PTNS_FUNCTION_MCMC_HEAT_BATH_H
#define TNBP_PTNS_FUNCTION_MCMC_HEAT_BATH_H
#include "tnbp/framework/graph_contraction.h"
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
    std::vector<IntT>& leaves,std::vector<int>& dimensions) {
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
    int failed=0;MPI_Allreduce(&bad,&failed,1,MPI_INT,MPI_MAX,comm);
    if(failed) throw std::invalid_argument("heat-bath graph/tensor/state/owner mismatch");
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
}


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

/** @brief Join two owner-local cached regions with the existing graph API.
 * Each external leg has a separate graph leaf, so common-neighbor bond fusion
 * cannot alter an independent cache. Shared legs are fused into one
 * composite bond on private operands; no numerical identity tensors are used.
 * The result stays on owner_a. Only copies of these two cached tensors and the
 * composite bond enter the contraction, never distant raw groups.
 * audit[0/1] counts real graph calls in environment construction/window evaluation;
 * audit[2] counts group preparation, [3] maximum region inputs per join,
 * [4] intermediate moves, [5] opposite-environment builds, [6] identity buffers,
 * [7] input clones in joins, [8] private boundary candidates, [9] binary graph
 * contractions. Counts are global logical operations. No hidden state.
 * Optional trace records {result owner, other owner, shared legs, dense MACs,
 * actual output elements}; it is populated only for explicit diagnostic calls.
 */
template<class TenT,class IntT>
void mcmc_contract_cached_tensors(context_handle_t<TenT>& ctx,
    const TenT* a,const std::vector<IntT>& axes_a,int owner_a,
    const TenT* b,const std::vector<IntT>& axes_b,int owner_b,
    const std::vector<tcapi::bond_dim_t<TenT>>& dimensions,MPI_Comm comm,
    std::vector<TenT>& result,std::vector<IntT>& axes_result,
    int phase=1,std::vector<std::size_t>* audit=nullptr,bool invalid_input=false,
    std::vector<std::vector<std::size_t>>* trace=nullptr) {
    int rank;MPI_Comm_rank(comm,&rank);
    int invalid=(invalid_input || (a!=nullptr)!=(rank==owner_a) || (b!=nullptr)!=(rank==owner_b)) ? 1:0;
    int failed=0;MPI_Allreduce(&invalid,&failed,1,MPI_INT,MPI_MAX,comm);
    if(failed) throw std::invalid_argument("cached region is not owner-local");
    std::set<IntT> shared;
    for(auto token:axes_a) if(std::find(axes_b.begin(),axes_b.end(),token)!=axes_b.end()) shared.insert(token);
    if(shared.empty()) throw std::invalid_argument("cached regions have no connecting edge");
    // Collapse all shared legs into one composite graph bond on private copies.
    // Both ends use the same token order; external cache axes stay independent.
    using Dim=tcapi::bond_dim_t<TenT>;
    Dim fused=1;
    for(auto token:shared) {
        const auto d=dimensions.at(static_cast<std::size_t>(token));
        if(d<1 || fused>std::numeric_limits<Dim>::max()/d)
            throw std::overflow_error("heat-bath composite bond dimension");
        fused*=d;
    }
    const IntT left=0,right=1;IntT next=2;
    std::map<IntT,IntT> port_token;
    std::vector<std::pair<IntT,IntT>> graph;
    std::vector<TenT> tensors;std::vector<IntT> ids;
    std::map<IntT,IntT> owners{{left,static_cast<IntT>(owner_a)},{right,static_cast<IntT>(owner_b)}};
    auto append=[&](const TenT* source,const std::vector<IntT>& axes,IntT id) {
        if(!source) return;
        tensors.push_back(tcapi::copy(ctx,*source));ids.push_back(id);
        auto remaining=axes;
        for(auto token:shared) if(token!=*shared.begin()) {
            const auto kept=std::find(remaining.begin(),remaining.end(),*shared.begin())-remaining.begin();
            const auto removed=std::find(remaining.begin(),remaining.end(),token)-remaining.begin();
            detail::fuse_graph_bond_axes(ctx,tensors.back(),kept,removed);
            remaining.erase(remaining.begin()+removed);
        }
    };
    append(a,axes_a,left);append(b,axes_b,right);
    std::set<IntT> tokens(axes_a.begin(),axes_a.end());tokens.insert(axes_b.begin(),axes_b.end());
    for(auto token:tokens) {
        if(shared.count(token)) {
            if(token==*shared.begin()) graph.emplace_back(left,right);
        } else {
            const auto port=next++;port_token[port]=token;
            graph.emplace_back(std::find(axes_a.begin(),axes_a.end(),token)!=axes_a.end() ? left:right,port);
        }
    }
    std::vector<IntT> internal{left,right};
    auto schedule=extract_induced_edges(graph,internal);
    if(audit) {
        audit->at(static_cast<std::size_t>(phase))++;audit->at(3)=std::max(audit->at(3),std::size_t(2));
        audit->at(7)+=2; // one input clone per immutable operand
        audit->at(9)+=internal.size()-1; // connected internal vertices contracted away
    }
    graph_tensor_network_contraction(ctx,graph,tensors,ids,owners,schedule,
        std::map<IntT,int>{{left,1}},std::map<IntT,IntT>{},comm);
    axes_result.clear();
    for(auto edge:graph) {
        if(edge.first!=left && edge.second!=left) throw std::runtime_error("cached join left an internal edge");
        axes_result.push_back(port_token.at(edge.first==left ? edge.second:edge.first));
    }
    result.clear();
    if(rank==owner_a) {
        auto index=std::find(ids.begin(),ids.end(),left)-ids.begin();
        result.push_back(std::move(tensors.at(index)));
    }
    mcmc_canonicalize_cache(ctx,result,axes_result);
    if(trace) {
        std::size_t work=1,elements=1;
        for(auto token:tokens) work*=static_cast<std::size_t>(dimensions.at(token));
        unsigned long long actual_elements=1;
        if(rank==owner_a) for(auto d:tcapi::shape(ctx,result.at(0))) actual_elements*=d;
        MPI_Bcast(&actual_elements,1,MPI_UNSIGNED_LONG_LONG,owner_a,comm);
        elements=static_cast<std::size_t>(actual_elements);
        trace->push_back({static_cast<std::size_t>(owner_a),static_cast<std::size_t>(owner_b),shared.size(),work,elements});
    }
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
    std::vector<std::vector<int>>* boundary_owners=nullptr) {
    int rank;MPI_Comm_rank(comm,&rank);
    std::vector<std::pair<IntT,IntT>> graph;std::vector<TenT> raw;
    std::vector<IntT> ids,leaves;std::map<IntT,IntT> ranks;
    mcmc_prepare_heat_bath_graph(ctx,edges,W,site_idx,owners,upper,upper_ids,lower,lower_ids,
        line,groups,states,comm,graph,raw,ids,ranks,leaves,physical_dimensions);
    std::vector<IntT> required;
    for(const auto& group:groups) required.insert(required.end(),group.begin(),group.end());
    detail::validate_graph_tensors(ctx,required,graph,raw,ids,ranks,comm);
    std::vector<unsigned long long> dimensions(graph.size(),0);
    for(std::size_t j=0;j<ids.size();++j) {
        auto bonds=GetSurroundingBondIndex(ids[j],graph);auto shape=tcapi::shape(ctx,raw[j]);
        for(std::size_t k=0;k<bonds.size();++k) dimensions[bonds[k]]=static_cast<unsigned long long>(shape[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE,dimensions.data(),static_cast<int>(dimensions.size()),MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
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
}

/** @brief Attach one boundary before the target, then the remaining boundary.
 * The pre-target intermediate is common to all physical states. On a square
 * lattice this is L B, (L B) T(s), ((L B) T(s)) U (or its mirror).
 * More general groups absorb additional connected members only when needed.
 * Target entry transfers the result to its owner; all later joins keep it there.
 * Optional steps records executed member indices, with -1 for target entry.
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
    int state_offset=0,std::vector<std::vector<std::size_t>>* trace=nullptr) {
    int rank;MPI_Comm_rank(comm,&rank);
    auto connected=[](const std::vector<IntT>& a,const std::vector<IntT>& b) {
        for(auto token:a) if(std::find(b.begin(),b.end(),token)!=b.end()) return true;
        return false;
    };
    std::vector<bool> used(boundary.size(),false);
    std::vector<TenT> common;auto common_axes=side_axes;
    const std::vector<TenT>* input=&side;bool have=has_side;int common_owner=side_owner;
    auto attach_boundary=[&](std::size_t k) {
        if(have) {
            std::vector<TenT> joined;std::vector<IntT> joined_axes;
            mcmc_contract_cached_regions(ctx,boundary[k],boundary_axes[k],boundary_owners[k],
                *input,common_axes,common_owner,dimensions,comm,joined,joined_axes,phase,audit,trace);
            common=std::move(joined);common_axes=std::move(joined_axes);input=&common;
        } else {
            // The raw boundary stays immutable and alive for every candidate.
            // The graph join makes its own private copy when it is used.
            input=&boundary[k];common_axes=boundary_axes[k];
        }
        common_owner=boundary_owners[k];have=true;used[k]=true;
        if(steps) steps->push_back(static_cast<int>(k));
    };
    // Select the first boundary connected to the advancing side. At an open
    // end start from a boundary adjacent to the target instead.
    for(std::size_t k=0;k<boundary.size();++k)
        if(connected(has_side ? side_axes:target_axes,boundary_axes[k])) {attach_boundary(k);break;}
    while(have && !connected(common_axes,target_axes)) {
        std::size_t k=0;for(;k<boundary.size();++k)
            if(!used[k] && connected(common_axes,boundary_axes[k])) break;
        if(k==boundary.size()) throw std::invalid_argument("group cannot reach target from environment");
        attach_boundary(k);
    }
    // Determine remaining connected absorption order once, independent of state.
    auto final_axes=target_axes;
    if(have) {
        for(auto token:common_axes) {
            auto it=std::find(final_axes.begin(),final_axes.end(),token);
            if(it==final_axes.end()) final_axes.push_back(token);else final_axes.erase(it);
        }
    }
    std::vector<std::size_t> remaining;
    for(std::size_t n=std::count(used.begin(),used.end(),false);n>0;--n) {
        std::size_t k=0;for(;k<boundary.size();++k)
            if(!used[k] && connected(final_axes,boundary_axes[k])) break;
        if(k==boundary.size()) throw std::invalid_argument("disconnected heat-bath group");
        used[k]=true;remaining.push_back(k);
        for(auto token:boundary_axes[k]) {
            auto it=std::find(final_axes.begin(),final_axes.end(),token);
            if(it==final_axes.end()) final_axes.push_back(token);else final_axes.erase(it);
        }
    }
    candidates.clear();candidates.resize(count);
    for(int state=0;state<count;++state) {
        if(have) mcmc_contract_cached_tensors(ctx,rank==target_owner ? &target.at(state+state_offset):nullptr,
            target_axes,target_owner,input->empty() ? nullptr:&input->at(0),common_axes,common_owner,
            dimensions,comm,candidates[state],candidate_axes,phase,audit,false,trace);
        else {
            if(rank==target_owner) candidates[state].push_back(tcapi::copy(ctx,target.at(state+state_offset)));
            candidate_axes=target_axes;if(audit && phase==1) audit->at(8)++;
        }
        if(steps) steps->push_back(-1);
        for(auto k:remaining) {
            std::vector<TenT> joined;std::vector<IntT> joined_axes;
            mcmc_contract_cached_regions(ctx,candidates[state],candidate_axes,target_owner,
                boundary[k],boundary_axes[k],boundary_owners[k],dimensions,comm,joined,joined_axes,phase,audit,trace);
            candidates[state]=std::move(joined);candidate_axes=std::move(joined_axes);
            if(steps) steps->push_back(static_cast<int>(k));
        }
    }
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
    const std::vector<std::vector<int>>* boundary_owners=nullptr) {
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
            dimensions,comm,candidate,suffix_axes[pos-1],0,audit,nullptr,states[j]);
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
    const std::vector<std::vector<IntT>>& second_boundary_axes={},const std::vector<int>& second_boundary_owners={}) {
    using Real=tcapi::real_t<TenT>;int rank;MPI_Comm_rank(comm,&rank);
    if(first_dim<1 || second_dim<1 || static_cast<std::size_t>(first_dim)>
       static_cast<std::size_t>(std::numeric_limits<int>::max())/static_cast<std::size_t>(second_dim))
        throw std::overflow_error("heat-bath candidate count");
    mcmc_prepare_group_candidates(ctx,first,first_axes,first_owner,first_dim,
        first_boundary,first_boundary_axes,first_boundary_owners,prefix,prefix_axes,has_prefix,prefix_owner,
        dimensions,comm,candidates,candidate_axes,1,audit);
    std::vector<std::vector<TenT>> right;std::vector<IntT> right_axes;
    if(pair) mcmc_prepare_group_candidates(ctx,second,second_axes,second_owner,second_dim,
        second_boundary,second_boundary_axes,second_boundary_owners,suffix,suffix_axes,has_suffix,suffix_owner,
        dimensions,comm,right,right_axes,1,audit);
    amplitudes.assign(static_cast<std::size_t>(first_dim)*second_dim,Real(0));int bad=0;
    for(int t=0;t<second_dim;++t) for(int s=0;s<first_dim;++s) {
        std::vector<TenT> evaluated;std::vector<IntT> evaluated_axes;
        if(pair || has_suffix) mcmc_contract_cached_regions(ctx,candidates[s],candidate_axes,first_owner,
            pair ? right[t]:suffix,pair ? right_axes:suffix_axes,pair ? second_owner:suffix_owner,
            dimensions,comm,evaluated,evaluated_axes,1,audit);
        else evaluated_axes=candidate_axes;
        if(!evaluated_axes.empty()) throw std::runtime_error("window evaluation did not close the cached boundary");
        if(rank==first_owner) try {
            const auto a=std::abs(tcapi::get_elem(ctx,(pair || has_suffix) ? evaluated.at(0):candidates[s].at(0),tcapi::elem_coors_t<TenT>{}));
            if(!std::isfinite(a)) bad=1;
            amplitudes[static_cast<std::size_t>(s)+static_cast<std::size_t>(first_dim)*t]=a;
        } catch(const std::exception&) {bad=1;}
    }
    MPI_Bcast(&bad,1,MPI_INT,first_owner,comm);
    if(bad) throw std::runtime_error("invalid cached heat-bath amplitude");
    MPI_Bcast(amplitudes.data(),static_cast<int>(amplitudes.size()),GetMpiType<Real>(),first_owner,comm);
    mcmc_heat_bath_weights(amplitudes,weights);
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

/** @brief Distributed full-state heat bath with both side environments reused.
 * Prepare raw members and proper opposite-side caches once. Each candidate contracts
 * only window groups and cached sides. Move the selected leading intermediate
 * forward by one group; the overlapping group remains explicit for the next
 * window. Reverse sweeps exchange the roles of left and right. No rank-zero
 * gather, no class, and no candidate-dependent exterior reconstruction.
 */
template<class TenT,class IntT,class RNG>
void mcmc_update_line_heat_bath(context_handle_t<TenT>& ctx,
    const std::vector<std::pair<IntT,IntT>>& edges,const std::vector<TenT>& W,
    const std::vector<IntT>& site_idx,const std::map<IntT,IntT>& owners,
    const std::vector<TenT>& upper,const std::vector<IntT>& upper_ids,
    std::vector<TenT>& target,const std::vector<IntT>& target_ids,
    const std::vector<TenT>& lower,const std::vector<IntT>& lower_ids,
    const std::vector<IntT>& line,const std::vector<std::vector<IntT>>& groups,
    std::vector<IntT>& states,MPI_Comm comm,RNG& rng,int direction,int width,
    std::vector<std::size_t>* audit=nullptr) {
    if(audit) audit->assign(10,0);
    if(width!=1 && width!=2) throw std::invalid_argument("update_sites must be 1 or 2");
    if(line.empty()) throw std::invalid_argument("empty heat-bath line");
    if(width==2 && line.size()==1) return;
    int rank;MPI_Comm_rank(comm,&rank);int bad=0,failed=0;
    for(auto site:line) if(owners.at(site)==rank) {
        const auto it=std::find(target_ids.begin(),target_ids.end(),site);
        if(it==target_ids.end() || static_cast<std::size_t>(it-target_ids.begin())>=target.size()) bad=1;
    }
    MPI_Allreduce(&bad,&failed,1,MPI_INT,MPI_MAX,comm);
    if(failed) throw std::invalid_argument("missing heat-bath commit tensor");
    std::vector<std::vector<TenT>> slices;std::vector<std::vector<IntT>> axes;
    std::vector<std::vector<std::vector<TenT>>> boundary;
    std::vector<std::vector<std::vector<IntT>>> boundary_axes;
    std::vector<std::vector<int>> boundary_owners;
    std::vector<int> dims;std::vector<tcapi::bond_dim_t<TenT>> bond_dims;
    mcmc_prepare_heat_bath_cache(ctx,edges,W,site_idx,owners,upper,upper_ids,lower,lower_ids,
        line,groups,states,comm,slices,axes,dims,bond_dims,audit,&boundary,&boundary_axes,&boundary_owners);
    std::vector<int> selected(line.size()),group_owners(line.size());
    for(std::size_t j=0;j<line.size();++j) {
        group_owners[j]=static_cast<int>(owners.at(line[j]));
        if(rank==group_owners[j]) selected[j]=static_cast<int>(states.at(std::find(site_idx.begin(),site_idx.end(),line[j])-site_idx.begin()));
        MPI_Bcast(&selected[j],1,MPI_INT,group_owners[j],comm);
    }
    std::vector<std::size_t> order(line.size());std::iota(order.begin(),order.end(),0);
    if(direction) std::reverse(order.begin(),order.end());
    std::vector<std::vector<TenT>> suffix;std::vector<std::vector<IntT>> suffix_axes;
    mcmc_build_heat_bath_environments(ctx,slices,axes,group_owners,selected,order,bond_dims,comm,suffix,suffix_axes,width,audit,&boundary,&boundary_axes,&boundary_owners);
    std::vector<TenT> prefix;std::vector<IntT> prefix_axes;
    const std::vector<TenT> empty;const std::vector<IntT> empty_axes;
    for(std::size_t pos=0;pos+static_cast<std::size_t>(width)<=line.size();++pos) {
        auto a=order[pos],b=width==2 ? order[pos+1]:a;const auto end=pos+static_cast<std::size_t>(width);
        std::vector<tcapi::real_t<TenT>> amplitudes,weights;
        std::vector<std::vector<TenT>> candidates;std::vector<IntT> candidate_axes;
        mcmc_evaluate_window_weights(ctx,slices[a],axes[a],group_owners[a],dims[a],
            width==2 ? slices[b]:empty,width==2 ? axes[b]:empty_axes,group_owners[b],width==2 ? dims[b]:1,width==2,
            prefix,prefix_axes,pos>0,pos>0 ? group_owners[order[pos-1]]:group_owners[a],
            suffix[end],suffix_axes[end],end<line.size(),end<line.size() ? group_owners[order[end]]:group_owners[a],
            bond_dims,comm,amplitudes,weights,candidates,candidate_axes,audit,
            boundary[a],boundary_axes[a],boundary_owners[a],boundary[b],boundary_axes[b],boundary_owners[b]);
        unsigned long long chosen=0;
        if(rank==group_owners[a]) chosen=static_cast<unsigned long long>(mcmc_draw_heat_bath(weights,rng));
        MPI_Bcast(&chosen,1,MPI_UNSIGNED_LONG_LONG,group_owners[a],comm);
        selected[a]=static_cast<int>(chosen%static_cast<unsigned long long>(dims[a]));
        if(width==2) selected[b]=static_cast<int>(chosen/static_cast<unsigned long long>(dims[a]));
        mcmc_absorb_group(static_cast<std::size_t>(selected[a]),candidates,candidate_axes,prefix,prefix_axes,audit);
    }
    mcmc_commit_window_state(slices,site_idx,owners,line,selected,rank,target,target_ids,states);
}
} // namespace tnbp
#endif
