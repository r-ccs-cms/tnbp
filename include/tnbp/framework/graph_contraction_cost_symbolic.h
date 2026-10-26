/// This file is a part of r-ccs-cms/tnbp
#ifndef TNBP_FRAMEWORK_GRAPH_CONTRACTION_COST_SYMBOLIC_H
#define TNBP_FRAMEWORK_GRAPH_CONTRACTION_COST_SYMBOLIC_H
#include "tnbp/framework/graph_contraction_cost.h"
#include <string>
#include <tuple>

namespace tnbp {
// An axis is D^a chi^b; coefficients belong to the cost polynomial.
using graph_symbolic_power = std::pair<unsigned,unsigned>;
using graph_symbolic_polynomial = std::map<graph_symbolic_power,std::uint64_t>;
template<class Site> using graph_symbolic_shapes = std::map<Site,std::vector<graph_symbolic_power>>;
template<class Site> using graph_string_shapes = std::map<Site,std::vector<std::string>>;
// total MACs, step MACs, step output sizes, final edges/shapes, shared graph operations.
template<class Site> using graph_symbolic_cost_result = std::tuple<
    graph_symbolic_polynomial,std::vector<graph_symbolic_power>,std::vector<graph_symbolic_power>,
    std::vector<std::pair<Site,Site>>,graph_symbolic_shapes<Site>,
    std::vector<detail::GraphContractionPlan<Site>>>;
namespace detail {
inline graph_symbolic_power symbolic_multiply(graph_symbolic_power a,graph_symbolic_power b) {
    if(b.first>std::numeric_limits<unsigned>::max()-a.first ||
       b.second>std::numeric_limits<unsigned>::max()-a.second)
        throw std::overflow_error("symbolic dimension exponent overflow");
    return {a.first+b.first,a.second+b.second};
}
inline graph_symbolic_power parse_symbolic_axis(const std::string& s) {
    if(s=="1") return {0,0};
    if(s.empty()) throw std::invalid_argument("empty symbolic dimension");
    graph_symbolic_power p{0,0};
    for(char c:s) {
        if(c=='D') p=symbolic_multiply(p,{1,0});
        else if(c=='x') p=symbolic_multiply(p,{0,1});
        else throw std::invalid_argument("symbolic dimension expects D, x, or 1");
    }
    return p;
}
template<class Site>
graph_symbolic_shapes<Site> parse_graph_symbolic_shapes(
    const std::vector<std::pair<Site,Site>>& edges,const graph_string_shapes<Site>& strings) {
    std::map<Site,std::size_t> degree,axis;
    for(auto edge:edges){++degree[edge.first];++degree[edge.second];}
    graph_symbolic_shapes<Site> shapes;
    for(const auto& entry:strings) {
        auto& dst=shapes[entry.first];
        for(const auto& s:entry.second) dst.push_back(parse_symbolic_axis(s));
        if(dst.size()!=degree[entry.first])
            throw std::invalid_argument("symbolic shape/incident edge count mismatch");
    }
    std::set<std::pair<Site,Site>> unique;
    for(const auto& [a,b]:edges) {
        if(a==b || !unique.insert(make_edge(a,b)).second || !shapes.count(a) || !shapes.count(b))
            throw std::invalid_argument("invalid symbolic graph");
        // Incident axes follow edge order; no repeated full-edge scans.
        if(shapes.at(a).at(axis[a]++)!=shapes.at(b).at(axis[b]++))
            throw std::invalid_argument("symbolic bond dimensions disagree");
    }
    return shapes;
}
template<class Index>
void symbolic_fuse(std::vector<graph_symbolic_power>& shape,std::vector<Index>& bonds,
                   std::size_t kept,std::size_t removed) {
    const auto a=std::find(bonds.begin(),bonds.end(),kept)-bonds.begin();
    const auto b=std::find(bonds.begin(),bonds.end(),removed)-bonds.begin();
    shape.at(a)=symbolic_multiply(shape.at(a),shape.at(b));
    shape.erase(shape.begin()+b); bonds.erase(bonds.begin()+b);
}
template<class Site>
std::pair<graph_symbolic_power,graph_symbolic_power> graph_symbolic_step(
    std::vector<std::pair<Site,Site>>& edges,graph_symbolic_shapes<Site>& shapes,
    Site a,Site b,Site c,const GraphContractionPlan<Site>* prepared=nullptr) {
    const auto generated=prepared ? GraphContractionPlan<Site>{}:make_graph_contraction_plan(edges,a,b,c);
    const auto& plan=prepared ? *prepared:generated;
    if(c!=a && shapes.count(c)) throw std::invalid_argument("symbolic result label collision");
    std::map<std::size_t,graph_symbolic_power> dimensions;
    for(std::size_t j=0;j<plan.bonds_a.size();++j) dimensions.emplace(plan.bonds_a[j],shapes.at(a).at(j));
    for(std::size_t j=0;j<plan.bonds_b.size();++j) dimensions.emplace(plan.bonds_b[j],shapes.at(b).at(j));
    graph_symbolic_power mac{0,0},elements{0,0};
    for(auto dim:dimensions) mac=symbolic_multiply(mac,dim.second);
    std::vector<graph_symbolic_power> out;
    for(auto e:plan.result_bonds) out.push_back(dimensions.at(e));
    auto bonds=plan.result_bonds;
    for(auto f:plan.fusions) {
        symbolic_fuse(out,bonds,f.kept_edge,f.removed_edge);
        auto neighbor_bonds=GetSurroundingBondIndex(f.neighbor,edges);
        symbolic_fuse(shapes.at(f.neighbor),neighbor_bonds,f.kept_edge,f.removed_edge);
    }
    for(auto dim:out) elements=symbolic_multiply(elements,dim);
    shapes.erase(a);shapes.erase(b);shapes.emplace(c,std::move(out));edges=plan.edges;
    return {mac,elements};
}
inline std::uint64_t symbolic_power_value(graph_symbolic_power p,std::uint64_t d,std::uint64_t chi) {
    if(!d || !chi) throw std::invalid_argument("model dimensions must be positive");
    auto pow=[](std::uint64_t x,unsigned n) {
        std::uint64_t v=1;
        while(n) {if(n&1) v=graph_cost_product(v,x);n>>=1;if(n)x=graph_cost_product(x,x);}
        return v;
    };
    return graph_cost_product(pow(d,p.first),pow(chi,p.second));
}
} // namespace detail

/** Dense cost with dimensions expressed as strings (D, x=chi, or 1).
 * Axes follow the original incident-edge order, exactly as the numerical cost
 * evaluator. The input graph/shapes are never changed, including on failure.
 * Steps are explicit {site_a,site_b,site_c}, not an edge-priority schedule.
 * All graph vertices, including common neighbors, require shapes.
 */
template<class Site>
graph_symbolic_cost_result<Site> graph_schedule_cost_symbolic(
    const std::vector<std::pair<Site,Site>>& edges,const graph_string_shapes<Site>& strings,
    const std::vector<std::array<Site,3>>& steps) {
    auto shapes=detail::parse_graph_symbolic_shapes(edges,strings);auto graph=edges;
    graph_symbolic_polynomial total;
    std::vector<graph_symbolic_power> macs,sizes;
    std::vector<detail::GraphContractionPlan<Site>> operations;
    for(auto step:steps) {
        operations.push_back(detail::make_graph_contraction_plan(graph,step[0],step[1],step[2]));
        const auto result=detail::graph_symbolic_step(graph,shapes,step[0],step[1],step[2],&operations.back());
        auto& coefficient=total[result.first];
        if(coefficient==std::numeric_limits<std::uint64_t>::max()) throw std::overflow_error("symbolic coefficient overflow");
        ++coefficient;macs.push_back(result.first);sizes.push_back(result.second);
    }
    return {std::move(total),std::move(macs),std::move(sizes),std::move(graph),std::move(shapes),std::move(operations)};
}
template<class Site>
graph_symbolic_cost_result<Site> graph_contraction_cost_symbolic(
    const std::vector<std::pair<Site,Site>>& edges,const graph_string_shapes<Site>& strings,
    Site a,Site b,Site c) {
    return graph_schedule_cost_symbolic(edges,strings,std::vector<std::array<Site,3>>{{a,b,c}});
}
inline std::uint64_t graph_symbolic_cost_value(const graph_symbolic_polynomial& polynomial,
                                             std::uint64_t d,std::uint64_t chi) {
    if(!d || !chi) throw std::invalid_argument("model dimensions must be positive");
    std::uint64_t value=0;
    for(auto term:polynomial) {
        const auto n=detail::graph_cost_product(term.second,detail::symbolic_power_value(term.first,d,chi));
        if(n>std::numeric_limits<std::uint64_t>::max()-value) throw std::overflow_error("symbolic cost value overflow");
        value+=n;
    }
    return value;
}
} // namespace tnbp
#endif
