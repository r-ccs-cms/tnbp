/// This file is a part of r-ccs-cms/tnbp
#ifndef TNBP_FRAMEWORK_GRAPH_SYMBOLIC_REGION_PLAN_H
#define TNBP_FRAMEWORK_GRAPH_SYMBOLIC_REGION_PLAN_H
#include "tnbp/framework/graph_contraction_cost_symbolic.h"
#include <functional>
namespace tnbp {
// Explicit steps, candidate-weighted MAC polynomial, step result sizes,
// all required original owner ranks, step multiplicities.
template<class Site> using graph_symbolic_region_plan=std::tuple<
    std::vector<std::array<Site,3>>,graph_symbolic_polynomial,
    std::vector<graph_symbolic_power>,std::vector<int>,std::vector<std::uint64_t>,
    std::vector<detail::GraphContractionPlan<Site>>>;

namespace detail {
// Internal entry: the caller parsed and validated the shared input once.
template<class Site>
graph_symbolic_region_plan<Site> evaluate_graph_symbolic_region_prepared(
    const std::vector<std::pair<Site,Site>>& edges,graph_symbolic_shapes<Site> shapes,
    const std::vector<std::array<Site,3>>& steps,const std::map<Site,int>& owners,
    Site target,std::uint64_t candidates) {
    if(!candidates) throw std::invalid_argument("zero physical candidate count");
    if(!shapes.count(target)) throw std::invalid_argument("missing repeated target");
    auto graph=edges;auto ranks=owners;std::set<Site> repeated{target};std::set<int> required;
    graph_symbolic_polynomial cost;std::vector<graph_symbolic_power> sizes;
    std::vector<std::uint64_t> repeats;
    std::vector<detail::GraphContractionPlan<Site>> operations;
    auto need=[&](Site site) {
        auto it=ranks.find(site);
        if(it==ranks.end() || it->second<0) throw std::invalid_argument("missing/invalid region owner");
        required.insert(it->second);
    };
    // Even a zero-step region must retain its owner's presence.
    need(target);
    for(auto st:steps) {
        auto plan=detail::make_graph_contraction_plan(graph,st[0],st[1],st[2]);
        operations.push_back(plan);
        need(st[0]);need(st[1]);for(auto f:plan.fusions)need(f.neighbor);
        const bool dependent=repeated.count(st[0]) || repeated.count(st[1]);
        const auto n=dependent ? candidates:1;
        const auto result=detail::graph_symbolic_step(graph,shapes,st[0],st[1],st[2],&plan);
        auto& coefficient=cost[result.first];
        if(n>std::numeric_limits<std::uint64_t>::max()-coefficient) throw std::overflow_error("weighted cost overflow");
        coefficient+=n;sizes.push_back(result.second);repeats.push_back(n);
        const auto rank=ranks.at(st[0]);ranks.erase(st[0]);ranks.erase(st[1]);ranks[st[2]]=rank;
        repeated.erase(st[0]);repeated.erase(st[1]);if(dependent)repeated.insert(st[2]);
    }
    return {steps,std::move(cost),std::move(sizes),std::vector<int>(required.begin(),required.end()),std::move(repeats),std::move(operations)};
}

} // namespace detail

template<class Site>
graph_symbolic_region_plan<Site> evaluate_graph_symbolic_region(
    const std::vector<std::pair<Site,Site>>& edges,const graph_string_shapes<Site>& strings,
    const std::vector<std::array<Site,3>>& steps,const std::map<Site,int>& owners,
    Site target,std::uint64_t candidates) {
    return detail::evaluate_graph_symbolic_region_prepared(edges,
        detail::parse_graph_symbolic_shapes(edges,strings),steps,owners,target,candidates);
}

/** All connected binary trees on <=5 region leaves, with target-owner retention.
 * Region leaves are a subset of a full original graph; outside common neighbors
 * remain in the graph and are included in the returned communication footprint.
 * The baseline is returned separately and is never removed by model selection.
 * No numerical shape or MPI operation is needed.
 */
namespace detail {
template<class Site>
std::pair<std::vector<graph_symbolic_region_plan<Site>>,graph_symbolic_region_plan<Site>>
build_graph_symbolic_region_plans_prepared(const std::vector<std::pair<Site,Site>>& edges,
    const graph_symbolic_shapes<Site>& shapes,const std::vector<Site>& leaves,
    const std::map<Site,int>& owners,Site target,std::uint64_t candidates,
    const std::vector<std::array<Site,3>>& baseline,std::size_t max_leaves=5) {
    auto validate_schedule=[&](const std::vector<std::array<Site,3>>& steps) {
        std::set<Site> live(leaves.begin(),leaves.end());
        if(live.empty() || live.size()!=leaves.size() || !live.count(target) || steps.size()!=live.size()-1)
            throw std::invalid_argument("invalid region leaves/schedule length");
        for(auto st:steps) {
            if(!live.count(st[0]) || !live.count(st[1]) || st[0]==st[1] || st[2]!=st[0] || st[1]==target)
                throw std::invalid_argument("region schedule must retain target and original A label");
            live.erase(st[1]);
        }
        if(*live.begin()!=target) throw std::invalid_argument("region final owner is not target");
    };
    validate_schedule(baseline);
    auto fallback=evaluate_graph_symbolic_region_prepared(edges,shapes,baseline,owners,target,candidates);
    std::vector<graph_symbolic_region_plan<Site>> result;
    if(leaves.size()>std::min<std::size_t>(max_leaves,5)) return {result,std::move(fallback)};
    std::vector<std::array<Site,3>> steps;
    std::function<void(const std::vector<std::pair<Site,Site>>&,std::vector<Site>)> visit;
    visit=[&](const auto& graph,std::vector<Site> live) {
        if(live.size()==1) {result.push_back(evaluate_graph_symbolic_region_prepared(edges,shapes,steps,owners,target,candidates));return;}
        std::sort(live.begin(),live.end());
        for(std::size_t i=0;i<live.size();++i)for(std::size_t j=i+1;j<live.size();++j) {
            Site a=live[i],b=live[j];if(b==target)std::swap(a,b);
            if(find_edge(graph,a,b)==graph.end())continue;
            auto next=graph;graph_tensor_contraction_for_edges(next,a,b,a);
            auto remaining=live;remaining.erase(std::find(remaining.begin(),remaining.end(),b));
            steps.push_back({a,b,a});visit(next,std::move(remaining));steps.pop_back();
        }
    };
    visit(edges,leaves);
    return {std::move(result),std::move(fallback)};
}

} // namespace detail

template<class Site>
std::pair<std::vector<graph_symbolic_region_plan<Site>>,graph_symbolic_region_plan<Site>>
build_graph_symbolic_region_plans(const std::vector<std::pair<Site,Site>>& edges,
    const graph_string_shapes<Site>& strings,const std::vector<Site>& leaves,
    const std::map<Site,int>& owners,Site target,std::uint64_t candidates,
    const std::vector<std::array<Site,3>>& baseline,std::size_t max_leaves=5) {
    const auto shapes=detail::parse_graph_symbolic_shapes(edges,strings);
    return detail::build_graph_symbolic_region_plans_prepared(edges,shapes,leaves,owners,target,candidates,baseline,max_leaves);
}

/** Select a prepared plan with known model dimensions (not measured bond sizes).
 * -1 means baseline. Incomparable/overflowing model evaluation uses baseline.
 * Lexicographic explicit steps break exact cost/output-size ties.
 */
template<class Site>
int select_graph_symbolic_region_plan(const std::vector<graph_symbolic_region_plan<Site>>& plans,
    const graph_symbolic_region_plan<Site>& baseline,std::uint64_t d,std::uint64_t chi) {
    auto score=[&](const auto& plan) {
        auto total=graph_symbolic_cost_value(std::get<1>(plan),d,chi);std::uint64_t peak=0;
        for(auto p:std::get<2>(plan))peak=std::max(peak,detail::symbolic_power_value(p,d,chi));
        return std::make_pair(total,peak);
    };
    if(!d || !chi) return -1;
    try {
        auto best=score(baseline);int chosen=-1;
        for(std::size_t i=0;i<plans.size();++i) {
            const auto value=score(plans[i]);
            const auto& current=chosen<0 ? baseline:plans[chosen];
            if(value<best || (value==best && std::get<0>(plans[i])<std::get<0>(current))) {
                best=value;chosen=static_cast<int>(i);
            }
        }
        return chosen;
    } catch(const std::overflow_error&) {return -1;}
}
} // namespace tnbp
#endif
