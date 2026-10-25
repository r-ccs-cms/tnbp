/// This file is a part of r-ccs-cms/tnbp
/** @file graph_contraction_cost.h
 * @brief Dense graph contraction costs from ordered edges and global shapes.
 */
#ifndef TNBP_FRAMEWORK_GRAPH_CONTRACTION_COST_H
#define TNBP_FRAMEWORK_GRAPH_CONTRACTION_COST_H
#include "tnbp/framework/graph_contraction_plan.h"
#include <cstdint>
#include <limits>
#include <type_traits>

namespace tnbp {
namespace detail {
// Optional reports must not participate in template deduction (allow nullptr).
template<class T> using graph_cost_output=typename std::common_type<T>::type;

inline std::uint64_t graph_cost_product(std::uint64_t a,std::uint64_t b) {
    if(b && a>std::numeric_limits<std::uint64_t>::max()/b)
        throw std::overflow_error("graph contraction cost: element/MAC count overflow");
    return a*b;
}

template<class DimT>
std::uint64_t graph_cost_dimension(DimT d) {
    static_assert(std::is_integral<DimT>::value && sizeof(DimT)<=sizeof(std::uint64_t),
                  "graph cost dimensions must be integral and at most 64 bits");
    if(d<1) throw std::invalid_argument("graph contraction cost: nonpositive dimension");
    return static_cast<std::uint64_t>(d);
}

template<class DimT>
std::uint64_t graph_cost_elements(const std::vector<DimT>& shape) {
    std::uint64_t n=1;
    for(auto d:shape) n=graph_cost_product(n,graph_cost_dimension(d));
    return n;
}

// All shapes are global metadata. Extra sites are permitted only as scalars.
template<class IntT,class DimT>
std::vector<DimT> graph_cost_dimensions(
    const std::vector<std::pair<IntT,IntT>>& edges,
    const std::map<IntT,std::vector<DimT>>& shapes) {
    if(edges.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("graph contraction cost: too many edges");
    std::set<std::pair<IntT,IntT>> unique;
    for(auto edge:edges) {
        if(edge.first==edge.second || !unique.insert(make_edge(edge.first,edge.second)).second)
            throw std::invalid_argument("graph contraction cost: invalid simple graph");
        if(!shapes.count(edge.first) || !shapes.count(edge.second))
            throw std::invalid_argument("graph contraction cost: missing site shape");
    }
    std::vector<DimT> dimensions(edges.size());std::vector<bool> seen(edges.size(),false);
    for(const auto& entry:shapes) {
        const auto bonds=GetSurroundingBondIndex(entry.first,edges);
        if(entry.second.size()!=bonds.size())
            throw std::invalid_argument("graph contraction cost: shape/incident-edge rank mismatch");
        for(std::size_t axis=0;axis<bonds.size();++axis) {
            graph_cost_dimension(entry.second[axis]);
            const auto e=static_cast<std::size_t>(bonds[axis]);
            if(seen[e] && dimensions[e]!=entry.second[axis])
                throw std::invalid_argument("graph contraction cost: bond dimensions disagree");
            seen[e]=true;dimensions[e]=entry.second[axis];
        }
    }
    return dimensions;
}

template<class DimT,class IndexT>
void graph_cost_fuse_axes(std::vector<DimT>& shape,std::vector<IndexT>& bonds,
                          std::size_t kept_edge,std::size_t removed_edge) {
    const auto kept=static_cast<std::size_t>(std::find(bonds.begin(),bonds.end(),kept_edge)-bonds.begin());
    const auto removed=static_cast<std::size_t>(std::find(bonds.begin(),bonds.end(),removed_edge)-bonds.begin());
    if(kept==bonds.size() || removed==bonds.size())
        throw std::logic_error("graph contraction cost: missing fused axis");
    const auto d=graph_cost_product(graph_cost_dimension(shape.at(kept)),graph_cost_dimension(shape.at(removed)));
    if(d>static_cast<std::uint64_t>(std::numeric_limits<DimT>::max()))
        throw std::overflow_error("graph contraction cost: fused dimension overflow");
    shape[kept]=static_cast<DimT>(d);shape.erase(shape.begin()+removed);
    bonds.erase(bonds.begin()+removed);
}
// Mutates private working metadata only; callers validate the complete input
// once and publish it only on success. Exceptions may leave this work partial.
template<class IntT,class DimT>
std::uint64_t graph_cost_step(
    std::vector<std::pair<IntT,IntT>>& edges,
    std::map<IntT,std::vector<DimT>>& shapes,IntT site_a,IntT site_b,IntT site_c) {
    auto plan=detail::make_graph_contraction_plan(edges,site_a,site_b,site_c);
    if(site_c!=site_a && shapes.count(site_c))
        throw std::invalid_argument("graph contraction cost: result label already exists");
    // Only operand axes are needed after the public entry-point validation.
    std::map<std::size_t,DimT> dimensions;
    for(std::size_t axis=0;axis<plan.bonds_a.size();++axis)
        dimensions.emplace(plan.bonds_a[axis],shapes.at(site_a).at(axis));
    for(std::size_t axis=0;axis<plan.bonds_b.size();++axis)
        dimensions.emplace(plan.bonds_b[axis],shapes.at(site_b).at(axis));
    std::uint64_t macs=1;
    for(const auto& entry:dimensions)
        macs=graph_cost_product(macs,graph_cost_dimension(entry.second));
    std::vector<DimT> result_shape;
    for(auto e:plan.result_bonds) result_shape.push_back(dimensions.at(e));
    auto result_bonds=plan.result_bonds;
    for(const auto& fusion:plan.fusions) {
        detail::graph_cost_fuse_axes(result_shape,result_bonds,fusion.kept_edge,fusion.removed_edge);
        auto bonds=GetSurroundingBondIndex(fusion.neighbor,edges);
        detail::graph_cost_fuse_axes(shapes.at(fusion.neighbor),bonds,fusion.kept_edge,fusion.removed_edge);
    }
    shapes.erase(site_a);shapes.erase(site_b);
    shapes.emplace(site_c,std::move(result_shape));
    edges.swap(plan.edges);
    return macs;
}
} // namespace detail

/** @brief Estimate one graph_tensor_contraction and update edges/shapes only.
 *
 * Shapes are global, keyed by site label, in current incident-edge order. All
 * graph vertices need shapes, including common neighbors; isolated scalars may
 * be present with empty shapes. Dimensions must be positive integral values.
 * Uses the numerical operation's graph plan, including the kept/removed axis
 * roles and updates to common neighbors. site_c is site_a or a fresh label.
 *
 * Returns dense scalar multiply-accumulate count: product of dimensions on the
 * union of the operands' incident edges (the contracted edge counted once).
 * No tensor values, backend, MPI or owner map are needed. This counts MACs, not
 * real/complex FLOPs, copies, transpose traffic, workspace or communication.
 * Scalar shape {} has one element. Overflow throws; inputs are unchanged on
 * validation or arithmetic failure.
 */
template<class IntT,class DimT>
std::uint64_t graph_tensor_contraction_cost(
    std::vector<std::pair<IntT,IntT>>& edges,
    std::map<IntT,std::vector<DimT>>& shapes,IntT site_a,IntT site_b,IntT site_c) {
    detail::graph_cost_dimensions(edges,shapes);
    auto work_edges=edges;auto work_shapes=shapes;
    const auto macs=detail::graph_cost_step(work_edges,work_shapes,site_a,site_b,site_c);
    edges.swap(work_edges);shapes.swap(work_shapes);
    return macs;
}

/** @brief Evaluate exactly a graph_tensor_network_contraction edge schedule.
 *
 * Priority, tie-breaking, label_map and schedule rewriting share the numerical
 * executor's helpers. Each step can change graph edges, operand shapes, and
 * common-neighbor shapes. The schedule may be partial or empty; no graph-size
 * limit is imposed beyond the executor's index limits. No order search occurs.
 *
 * Returns total dense MACs and replaces edges/shapes with the final metadata.
 * Optional outputs have one entry per executed binary contraction (not per
 * original scheduled edge, since redundant scheduled edges are fused away):
 * - step_macs: dense MACs per contraction.
 * - step_elements: resulting tensor elements, after common-neighbor fusion.
 *   Their maximum is the largest result tensor, NOT peak allocated memory;
 *   input tensors, private copies and backend workspace are not included.
 * - steps: executed {site_a,site_b,site_c}.
 * - shape_history/edge_history: post-step snapshots, no initial entry.
 *
 * Inputs and optional outputs are unchanged if evaluation fails. Optional
 * output containers must be distinct. For comparing orders, pass copies of
 * the same initial edges/shapes to each evaluation.
 */
template<class IntT,class DimT,class LevelT=int>
std::uint64_t graph_tensor_network_contraction_cost(
    std::vector<std::pair<IntT,IntT>>& edges,
    std::map<IntT,std::vector<DimT>>& shapes,
    const std::vector<std::pair<IntT,IntT>>& edges_schedule,
    const std::map<IntT,LevelT>& site_priority={},const std::map<IntT,IntT>& label_map={},
    std::vector<std::uint64_t>* step_macs=nullptr,std::vector<std::uint64_t>* step_elements=nullptr,
    detail::graph_cost_output<std::vector<std::array<IntT,3>>>* steps=nullptr,
    detail::graph_cost_output<std::vector<std::map<IntT,std::vector<DimT>>>>* shape_history=nullptr,
    detail::graph_cost_output<std::vector<std::vector<std::pair<IntT,IntT>>>>* edge_history=nullptr) {
    if(step_macs && step_macs==step_elements)
        throw std::invalid_argument("graph contraction cost: aliased output containers");
    detail::graph_cost_dimensions(edges,shapes);
    auto work_edges=edges;auto work_shapes=shapes;auto schedule=edges_schedule;
    const auto priority=detail::graph_contraction_priorities(schedule,site_priority,label_map);
    std::vector<std::uint64_t> macs,elements;
    std::vector<std::array<IntT,3>> executed;
    std::vector<std::map<IntT,std::vector<DimT>>> shapes_trace;
    std::vector<std::vector<std::pair<IntT,IntT>>> edges_trace;
    std::uint64_t total=0;
    while(!schedule.empty()) {
        const auto op=detail::graph_contraction_sites(schedule.front(),priority,label_map);
        const auto cost=detail::graph_cost_step(work_edges,work_shapes,op[0],op[1],op[2]);
        if(total>std::numeric_limits<std::uint64_t>::max()-cost)
            throw std::overflow_error("graph contraction cost: total MAC count overflow");
        total+=cost;
        graph_tensor_contraction_for_edges(schedule,op[0],op[1],op[2]);
        if(step_macs) macs.push_back(cost);
        if(step_elements) elements.push_back(detail::graph_cost_elements(work_shapes.at(op[2])));
        if(steps) executed.push_back(op);
        if(shape_history) shapes_trace.push_back(work_shapes);
        if(edge_history) edges_trace.push_back(work_edges);
    }
    edges.swap(work_edges);shapes.swap(work_shapes);
    if(step_macs) step_macs->swap(macs);
    if(step_elements) step_elements->swap(elements);
    if(steps) steps->swap(executed);
    if(shape_history) shape_history->swap(shapes_trace);
    if(edge_history) edge_history->swap(edges_trace);
    return total;
}
} // namespace tnbp
#endif
