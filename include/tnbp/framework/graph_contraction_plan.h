/// This file is a part of r-ccs-cms/tnbp
/** @file graph_contraction_plan.h
 * @brief Shared graph updates and schedule selection, independent of tensors.
 */
#ifndef TNBP_FRAMEWORK_GRAPH_CONTRACTION_PLAN_H
#define TNBP_FRAMEWORK_GRAPH_CONTRACTION_PLAN_H
#include "tnbp/framework/graph.h"
#include <array>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
namespace tnbp {
namespace detail {
// Immutable graph plan. Edge indices refer to the original edge vector;
// fusions are ordered by the removed operand's incident-edge order.
template<class IntT>
struct GraphContractionPlan {
    struct Fusion { IntT neighbor; std::size_t kept_edge, removed_edge; };
    std::vector<std::size_t> bonds_a, bonds_b, result_bonds;
    std::vector<Fusion> fusions;
    std::vector<std::pair<IntT,IntT>> edges;
};

template<class IntT>
GraphContractionPlan<IntT> make_graph_contraction_plan(
    const std::vector<std::pair<IntT,IntT>>& edges, IntT a, IntT b, IntT c) {
    if(a==b || c==b || find_edge(edges,a,b)==edges.end())
        throw std::invalid_argument("graph contraction: invalid operands or result label");
    GraphContractionPlan<IntT> plan;
    std::map<IntT,std::size_t> neighbors_a;
    std::set<std::pair<IntT,IntT>> unique;
    for(std::size_t e=0;e<edges.size();++e) {
        const auto [x,y]=edges[e];
        if(x==y || !unique.insert(make_edge(x,y)).second ||
           (c!=a && (x==c || y==c)))
            throw std::invalid_argument("graph contraction: invalid graph or result collision");
        const bool at_a=x==a || y==a, at_b=x==b || y==b;
        if(at_a) {plan.bonds_a.push_back(e);neighbors_a.emplace(x==a ? y : x,e);}
        if(at_b) plan.bonds_b.push_back(e);
        if(at_a!=at_b) plan.result_bonds.push_back(e);
    }
    std::set<std::size_t> removed;
    for(auto e:plan.bonds_b) {
        const auto [x,y]=edges[e];const auto other=x==b ? y : x;
        if(other==a) {removed.insert(e);continue;}
        const auto kept=neighbors_a.find(other);
        if(kept!=neighbors_a.end()) {
            plan.fusions.push_back({other,kept->second,e});removed.insert(e);
        }
    }
    // Keep original vector positions relative to surviving edges. Shared bonds
    // keep A's edge; B-only bonds are relabeled in place, never appended/sorted.
    for(std::size_t e=0;e<edges.size();++e) {
        if(removed.count(e)) continue;
        auto [x,y]=edges[e];
        if(x==a || x==b || y==a || y==b) {
            if(x==a || x==b) x=c;
            if(y==a || y==b) y=c;
            plan.edges.push_back(make_edge(x,y));
        } else plan.edges.push_back(edges[e]);
    }
    return plan;
}

// Keep scheduling decisions shared by numerical execution and shape-only cost.
template<class IntT, class LevelT>
std::map<IntT,LevelT> graph_contraction_priorities(
    const std::vector<std::pair<IntT,IntT>>& schedule,
    const std::map<IntT,LevelT>& priority, const std::map<IntT,IntT>& labels) {
    auto result=priority;
    for(auto site:GetSiteIndexFromBond(schedule)) {
        if(result.find(site)!=result.end()) {
            if(labels.find(site)!=labels.end()) result[labels.at(site)]=result.at(site);
        } else {
            result[site]=0;
            if(labels.find(site)!=labels.end()) result[labels.at(site)]=0;
        }
    }
    return result;
}

template<class IntT, class LevelT>
std::array<IntT,3> graph_contraction_sites(const std::pair<IntT,IntT>& edge,
    const std::map<IntT,LevelT>& priority, const std::map<IntT,IntT>& labels) {
    const auto a=priority.at(edge.second)<=priority.at(edge.first) ? edge.first:edge.second;
    const auto b=a==edge.first ? edge.second:edge.first;
    return {a,b,labels.find(a)!=labels.end() ? labels.at(a):a};
}
} // namespace detail

/** @brief Apply exactly the graph update of graph_tensor_contraction. */
template<class IntT>
void graph_tensor_contraction_for_edges(std::vector<std::pair<IntT,IntT>>& edges,
    IntT site_a, IntT site_b, IntT site_c) {
    edges=detail::make_graph_contraction_plan(edges,site_a,site_b,site_c).edges;
}
} // namespace tnbp
#endif
