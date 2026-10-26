/// This file is a part of r-ccs-cms/tnbp
#ifndef TNBP_FRAMEWORK_GRAPH_REGION_CONTRACTION_PLAN_H
#define TNBP_FRAMEWORK_GRAPH_REGION_CONTRACTION_PLAN_H
#include "tnbp/framework/graph_contraction_cost.h"
#include <tuple>
#include <iterator>
#include <numeric>

namespace tnbp {
namespace detail {
template<class TokenT,class DimT>
void validate_region_axes(const std::vector<TokenT>& axes,const std::vector<DimT>& dims) {
    if(!std::is_sorted(axes.begin(),axes.end()) ||
       std::adjacent_find(axes.begin(),axes.end())!=axes.end())
        throw std::invalid_argument("region axes must be sorted and unique");
    for(auto t:axes) {
        if(t<0 || static_cast<std::size_t>(t)>=dims.size())
            throw std::invalid_argument("invalid region token");
        graph_cost_dimension(dims.at(static_cast<std::size_t>(t)));
    }
}
template<class TokenT>
bool regions_connected(const std::vector<TokenT>& a,const std::vector<TokenT>& b) {
    auto i=a.begin(),j=b.begin();
    while(i!=a.end() && j!=b.end()) {
        if(*i==*j) return true;
        if(*i<*j) ++i;else ++j;
    }
    return false;
}
inline std::uint64_t region_cost_sum(std::uint64_t a,std::uint64_t b) {
    if(a>std::numeric_limits<std::uint64_t>::max()-b)
        throw std::overflow_error("region total MAC overflow");
    return a+b;
}
} // namespace detail

/** Prepare exactly the small graph used to join two immutable cached regions.
 * Each external token gets its own metadata leaf. Shared tokens are fused in
 * ascending token order on both operands; fusion indices refer to evolving axes.
 * All outputs are distinct and unchanged on failure. No tensors or MPI needed.
 */
template<class TokenT,class DimT>
void make_region_join_metadata(const std::vector<TokenT>& a,const std::vector<TokenT>& b,
    const std::vector<DimT>& dims,std::vector<std::pair<int,int>>& edges,
    std::map<int,std::vector<DimT>>& shapes,
    std::vector<std::pair<std::size_t,std::size_t>>& fuse_a,
    std::vector<std::pair<std::size_t,std::size_t>>& fuse_b,std::vector<TokenT>& result_axes) {
    detail::validate_region_axes(a,dims);detail::validate_region_axes(b,dims);
    std::vector<TokenT> shared,tokens,out;
    std::set_intersection(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(shared));
    if(shared.empty()) throw std::invalid_argument("region join requires a shared bond");
    std::set_union(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(tokens));
    if(tokens.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()-2))
        throw std::overflow_error("too many region ports");
    std::uint64_t fused=1;
    for(auto t:shared) fused=detail::graph_cost_product(fused,detail::graph_cost_dimension(dims.at(t)));
    if(fused>static_cast<std::uint64_t>(std::numeric_limits<DimT>::max()))
        throw std::overflow_error("region composite dimension overflow");
    auto fusion=[&](std::vector<TokenT> axes) {
        std::vector<std::pair<std::size_t,std::size_t>> f;
        for(std::size_t k=1;k<shared.size();++k) {
            const auto kept=std::find(axes.begin(),axes.end(),shared[0])-axes.begin();
            const auto removed=std::find(axes.begin(),axes.end(),shared[k])-axes.begin();
            f.emplace_back(kept,removed);axes.erase(axes.begin()+removed);
        }
        return f;
    };
    auto fa=fusion(a),fb=fusion(b);
    std::vector<std::pair<int,int>> graph;std::map<int,std::vector<DimT>> sh;
    sh[0]={};sh[1]={};int next=2;
    for(auto t:tokens) {
        if(std::binary_search(shared.begin(),shared.end(),t)) {
            if(t!=shared[0]) continue;
            graph.emplace_back(0,1);sh[0].push_back(static_cast<DimT>(fused));sh[1].push_back(static_cast<DimT>(fused));
        } else {
            const int id=std::binary_search(a.begin(),a.end(),t)?0:1;
            graph.emplace_back(id,next);sh[id].push_back(dims.at(t));sh[next++]={dims.at(t)};out.push_back(t);
        }
    }
    edges.swap(graph);shapes.swap(sh);fuse_a.swap(fa);fuse_b.swap(fb);result_axes.swap(out);
}

namespace detail {
// Structural normalization is independent of dimension arithmetic, so an
// overflowing baseline does not prevent exploring another contraction tree.
inline std::vector<std::array<std::size_t,3>> normalize_region_plan(std::size_t n,std::size_t repeated,
    const std::vector<std::array<std::size_t,3>>& steps,bool orient) {
    if(!n || n>std::numeric_limits<std::size_t>::max()/2 || repeated>=n || steps.size()!=n-1)
        throw std::invalid_argument("invalid region plan size");
    const auto nodes=2*n-1;
    std::vector<std::vector<std::size_t>> leaves(nodes);
    std::vector<unsigned char> used(nodes,0),dep(nodes,0);
    for(std::size_t i=0;i<n;++i) {leaves[i]={i};dep[i]=i==repeated;}
    auto ops=steps;
    for(std::size_t k=0;k<steps.size();++k) {
        auto& op=ops[k];auto a=op[0],b=op[1];const auto c=n+k;
        if(op[2]!=c || a>=c || b>=c || a==b || used[a] || used[b])
            throw std::invalid_argument("region plan is not a binary tree");
        used[a]=used[b]=1;dep[c]=dep[a] || dep[b];
        std::set_union(leaves[a].begin(),leaves[a].end(),leaves[b].begin(),leaves[b].end(),std::back_inserter(leaves[c]));
        if(orient) {
            if(dep[b] || (!dep[a] && leaves[b].front()<leaves[a].front())) std::swap(op[0],op[1]);
        } else if(dep[b]) throw std::invalid_argument("repeated operand must retain its owner");
    }
    for(std::size_t i=0;i+1<nodes;++i) if(!used[i]) throw std::invalid_argument("unused region node");
    std::vector<std::size_t> remap(nodes,std::numeric_limits<std::size_t>::max());
    for(std::size_t i=0;i<n;++i) remap[i]=i;
    std::vector<std::array<std::size_t,3>> out;
    while(out.size()<steps.size()) {
        std::size_t best=nodes;
        for(std::size_t k=0;k<ops.size();++k) {
            const auto op=ops[k];
            if(remap[op[2]]!=std::numeric_limits<std::size_t>::max() ||
               remap[op[0]]==std::numeric_limits<std::size_t>::max() ||
               remap[op[1]]==std::numeric_limits<std::size_t>::max()) continue;
            if(best==nodes || std::tie(dep[op[2]],leaves[op[2]])<std::tie(dep[best],leaves[best])) best=op[2];
        }
        if(best==nodes) throw std::logic_error("region plan ordering failed");
        const auto op=ops[best-n];const auto c=n+out.size();remap[best]=c;
        out.push_back({remap[op[0]],remap[op[1]],c});
    }
    return out;
}
template<class TokenT,class DimT>
void validate_region_inputs(const std::vector<std::vector<TokenT>>& axes,const std::vector<DimT>& dims,
    const std::vector<int>& owners,std::size_t repeated) {
    if(axes.empty() || axes.size()!=owners.size() || repeated>=axes.size())
        throw std::invalid_argument("invalid region leaves");
    std::map<TokenT,int> occurrences;
    for(std::size_t i=0;i<axes.size();++i) {
        validate_region_axes(axes[i],dims);
        if(owners[i]<0) throw std::invalid_argument("negative region owner");
        for(auto t:axes[i]) if(++occurrences[t]>2) throw std::invalid_argument("region hyperedge");
    }
}
template<class TokenT,class DimT>
std::uint64_t region_join_cost(const std::vector<TokenT>& a,const std::vector<TokenT>& b,
    const std::vector<DimT>& dims,std::vector<TokenT>& out,std::uint64_t& elements) {
    std::vector<std::pair<int,int>> graph;std::map<int,std::vector<DimT>> shapes;
    std::vector<std::pair<std::size_t,std::size_t>> fa,fb;
    make_region_join_metadata(a,b,dims,graph,shapes,fa,fb,out);
    const auto cost=graph_tensor_contraction_cost(graph,shapes,0,1,0);
    elements=graph_cost_elements(shapes.at(0));return cost;
}
} // namespace detail

/** Describe a normalized plan. Outputs must be distinct; failures are atomic. */
template<class TokenT,class DimT>
void describe_region_contraction_plan(const std::vector<std::vector<TokenT>>& leaf_axes,
    const std::vector<DimT>& dims,const std::vector<int>& leaf_owners,std::size_t repeated,
    const std::vector<std::array<std::size_t,3>>& steps,
    std::vector<std::vector<TokenT>>& node_axes,std::vector<int>& node_owners,
    std::vector<unsigned char>& node_repeated,std::vector<std::size_t>& parent,
    std::vector<std::size_t>& common_frontier,std::vector<std::uint64_t>& step_macs,
    std::vector<std::uint64_t>& step_elements) {
    detail::validate_region_inputs(leaf_axes,dims,leaf_owners,repeated);
    const auto n=leaf_axes.size();
    if(detail::normalize_region_plan(n,repeated,steps,n<=5)!=steps)
        throw std::invalid_argument("region plan must be normalized");
    auto axes=leaf_axes;auto owners=leaf_owners;
    axes.resize(2*n-1);owners.resize(2*n-1);
    std::vector<unsigned char> dep(2*n-1,0);dep[repeated]=1;
    std::vector<std::size_t> par(2*n-1,std::numeric_limits<std::size_t>::max()),frontier;
    std::vector<std::uint64_t> macs,elements;
    for(auto op:steps) {
        std::uint64_t size=0;
        macs.push_back(detail::region_join_cost(axes[op[0]],axes[op[1]],dims,axes[op[2]],size));
        elements.push_back(size);owners[op[2]]=owners[op[0]];
        dep[op[2]]=dep[op[0]] || dep[op[1]];par[op[0]]=par[op[1]]=op[2];
    }
    for(std::size_t i=0;i+1<par.size();++i) if(!dep[i] && dep[par[i]]) frontier.push_back(i);
    node_axes.swap(axes);node_owners.swap(owners);node_repeated.swap(dep);parent.swap(par);
    common_frontier.swap(frontier);step_macs.swap(macs);step_elements.swap(elements);
}

/** Exact connected binary-tree search for <=5 leaves, otherwise validated
 * baseline fallback. Costs include reuse of all nodes independent of repeated.
 * The baseline topology is retained but its orientation is normalized for <=5.
 * Returns total scalar MACs. No communication, tensor values or new classes.
 */
template<class TokenT,class DimT>
std::uint64_t make_small_region_contraction_plan(const std::vector<std::vector<TokenT>>& leaf_axes,
    const std::vector<DimT>& dims,const std::vector<int>& owners,std::size_t repeated,
    std::uint64_t repetitions,const std::vector<std::array<std::size_t,3>>& fallback,
    std::vector<std::array<std::size_t,3>>& steps,std::uint64_t& max_elements,bool& used_fallback) {
    detail::validate_region_inputs(leaf_axes,dims,owners,repeated);
    if(!repetitions) throw std::invalid_argument("zero region repetitions");
    const auto n=leaf_axes.size();auto baseline=detail::normalize_region_plan(n,repeated,fallback,n<=5);
    // Validate baseline connectivity even if dimension arithmetic overflows.
    auto baseline_axes=leaf_axes;baseline_axes.resize(2*n-1);
    for(auto op:baseline) {
        if(!detail::regions_connected(baseline_axes[op[0]],baseline_axes[op[1]]))
            throw std::invalid_argument("disconnected baseline");
        std::set_symmetric_difference(baseline_axes[op[0]].begin(),baseline_axes[op[0]].end(),
            baseline_axes[op[1]].begin(),baseline_axes[op[1]].end(),std::back_inserter(baseline_axes[op[2]]));
    }
    std::vector<std::array<std::size_t,3>> best;
    std::uint64_t best_cost=0,best_size=0;bool found=false;
    auto consider=[&](const auto& ops,std::uint64_t total,std::uint64_t peak) {
        auto canonical=detail::normalize_region_plan(n,repeated,ops,n<=5);
        if(!found || std::tie(total,peak,canonical)<std::tie(best_cost,best_size,best)) {
            best=std::move(canonical);best_cost=total;best_size=peak;found=true;
        }
    };
    try {
        std::vector<std::vector<TokenT>> a;std::vector<int> o;std::vector<unsigned char> d;
        std::vector<std::size_t> par,front;std::vector<std::uint64_t> macs,sizes;
        describe_region_contraction_plan(leaf_axes,dims,owners,repeated,baseline,a,o,d,par,front,macs,sizes);
        std::uint64_t total=0,peak=0;
        for(std::size_t k=0;k<macs.size();++k) {
            total=detail::region_cost_sum(total,detail::graph_cost_product(macs[k],d[n+k]?repetitions:1));
            peak=std::max(peak,sizes[k]);
        }
        consider(baseline,total,peak);
    } catch(const std::overflow_error&) {if(n>5) throw;}
    if(n<=5) {
        auto axes=leaf_axes;axes.resize(2*n-1);
        std::vector<unsigned char> dep(2*n-1,0);dep[repeated]=1;
        std::vector<std::size_t> first(2*n-1);std::iota(first.begin(),first.begin()+n,0);
        std::vector<std::size_t> active(n);std::iota(active.begin(),active.end(),0);
        std::vector<std::array<std::size_t,3>> ops;
        auto search=[&](auto&& self,const std::vector<std::size_t>& live,std::uint64_t total,std::uint64_t peak)->void {
            if(live.size()==1) {consider(ops,total,peak);return;}
            for(std::size_t i=0;i<live.size();++i) for(std::size_t j=i+1;j<live.size();++j) {
                auto a=live[i],b=live[j];const auto c=n+ops.size();
                if(!detail::regions_connected(axes[a],axes[b])) continue;
                if(dep[b] || (!dep[a] && first[b]<first[a])) std::swap(a,b);
                std::uint64_t size=0,cost=0,next_total=0;
                try {
                    cost=detail::region_join_cost(axes[a],axes[b],dims,axes[c],size);
                    next_total=detail::region_cost_sum(total,detail::graph_cost_product(cost,(dep[a]||dep[b])?repetitions:1));
                } catch(const std::overflow_error&) {continue;}
                dep[c]=dep[a] || dep[b];first[c]=std::min(first[a],first[b]);
                std::vector<std::size_t> next;
                for(std::size_t k=0;k<live.size();++k) if(k!=i && k!=j) next.push_back(live[k]);
                next.push_back(c);ops.push_back({a,b,c});self(self,next,next_total,std::max(peak,size));ops.pop_back();
            }
        };
        search(search,active,0,0);
    }
    if(!found) throw std::overflow_error("all region contraction plans overflow");
    steps.swap(best);max_elements=best_size;used_fallback=n>5;return best_cost;
}
} // namespace tnbp
#endif
