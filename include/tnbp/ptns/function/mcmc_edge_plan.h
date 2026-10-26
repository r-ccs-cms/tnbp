/// This file is a part of r-ccs-cms/tnbp
#ifndef TNBP_PTNS_FUNCTION_MCMC_EDGE_PLAN_H
#define TNBP_PTNS_FUNCTION_MCMC_EDGE_PLAN_H
#include "tnbp/framework/graph_symbolic_region_plan.h"
#include <numeric>
namespace tnbp {
// The ports are the ordered factors of an edge after native common-neighbor
// fusion. They refer to edges in the original three-line graph, not fake leaves.
using mcmc_edge_ports=std::vector<std::vector<std::size_t>>;
template<class S> using mcmc_edge_step=std::tuple<std::vector<std::pair<S,S>>,
    mcmc_edge_ports,detail::GraphContractionPlan<S>,mcmc_edge_ports>;
// initial graph/ports, steps, participants, target, compiled step metadata,
// common-prefix length, initial sources, sources of the repeated suffix.
template<class S> using mcmc_edge_region_plan=std::tuple<
    std::vector<std::pair<S,S>>,mcmc_edge_ports,std::vector<std::array<S,3>>,std::vector<int>,S,
    std::vector<mcmc_edge_step<S>>,std::size_t,std::vector<S>,std::vector<S>>;
template<class S> using mcmc_edge_window_plan=std::tuple<
    mcmc_edge_region_plan<S>,mcmc_edge_region_plan<S>,mcmc_edge_region_plan<S>,std::vector<int>,std::array<int,2>>;
template<class S> using mcmc_edge_line_plan=std::tuple<std::vector<std::size_t>,
    std::vector<mcmc_edge_region_plan<S>>,std::vector<mcmc_edge_window_plan<S>>,std::vector<std::vector<int>>>;

template<class S>
std::pair<std::vector<std::pair<S,S>>,mcmc_edge_ports> mcmc_edge_quotient(
    const std::vector<std::pair<S,S>>& original,const std::map<S,S>& representatives) {
    std::vector<std::pair<S,S>> edges;mcmc_edge_ports ports;std::map<std::pair<S,S>,std::size_t> index;
    for(std::size_t e=0;e<original.size();++e) {
        const auto a=representatives.at(original[e].first),b=representatives.at(original[e].second);
        if(a==b)continue;
        const auto edge=make_edge(a,b);auto found=index.find(edge);
        if(found==index.end()) {index[edge]=edges.size();edges.push_back(edge);ports.push_back({e});}
        else ports[found->second].push_back(e);
    }
    return {std::move(edges),std::move(ports)};
}
template<class S>
mcmc_edge_ports mcmc_edge_updated_ports(const std::vector<std::pair<S,S>>& edges,
    const mcmc_edge_ports& ports,const detail::GraphContractionPlan<S>& plan,S a,S b) {
    auto changed=ports;std::set<std::size_t> removed;
    removed.insert(static_cast<std::size_t>(find_edge(edges,a,b)-edges.begin()));
    for(auto f:plan.fusions) {
        auto& out=changed.at(f.kept_edge);const auto& in=changed.at(f.removed_edge);
        out.insert(out.end(),in.begin(),in.end());removed.insert(f.removed_edge);
    }
    mcmc_edge_ports result;
    for(std::size_t e=0;e<ports.size();++e)if(!removed.count(e))result.push_back(std::move(changed[e]));
    return result;
}
template<class S>
std::vector<std::size_t> mcmc_edge_cut(const std::vector<std::pair<S,S>>& edges,
    const mcmc_edge_ports& ports,S site) {
    std::vector<std::size_t> result;
    for(auto e:GetSurroundingBondIndex(site,edges))result.insert(result.end(),ports.at(e).begin(),ports.at(e).end());
    std::sort(result.begin(),result.end());return result;
}

template<class S>
void mcmc_compile_edge_region(mcmc_edge_region_plan<S>& region) {
    auto graph=std::get<0>(region);auto ports=std::get<1>(region);
    const auto& steps=std::get<2>(region);const auto target=std::get<4>(region);
    auto& compiled=std::get<5>(region);compiled.clear();
    std::size_t common=0;bool repeated=false;
    for(auto st:steps) {
        repeated=repeated || st[0]==target || st[1]==target;
        if(!repeated)++common;
        auto plan=detail::make_graph_contraction_plan(graph,st[0],st[1],st[2]);
        auto next=mcmc_edge_updated_ports(graph,ports,plan,st[0],st[1]);
        compiled.emplace_back(graph,ports,plan,next);graph=plan.edges;ports=std::move(next);
    }
    std::get<6>(region)=common;
    auto sources=[&](std::size_t start) {
        std::set<S> required,produced;
        for(std::size_t k=start;k<steps.size();++k) {
            const auto st=steps[k];std::vector<S> used{st[0],st[1]};
            for(auto f:std::get<2>(compiled[k]).fusions)used.push_back(f.neighbor);
            for(auto site:used)if(!produced.count(site))required.insert(site);
            produced.insert(st[2]);
        }
        if(start==steps.size())required.insert(target);
        return std::vector<S>(required.begin(),required.end());
    };
    std::get<7>(region)=sources(0);std::get<8>(region)=sources(common);
}

template<class S>
std::vector<std::array<S,3>> mcmc_edge_baseline(std::vector<std::pair<S,S>> graph,
    const std::vector<S>& leaves,S target,bool has_side,S side) {
    std::vector<std::array<S,3>> ops;std::set<S> remaining(leaves.begin(),leaves.end());
    remaining.erase(target);S current=target;
    auto join=[&](S a,S b) {ops.push_back({a,b,a});graph_tensor_contraction_for_edges(graph,a,b,a);remaining.erase(b);};
    bool seeded=false;
    for(auto site:leaves)if(site!=target && (!has_side || site!=side) &&
        find_edge(graph,has_side?side:target,site)!=graph.end()) {
        if(has_side) {current=site;join(site,side);}else current=site;
        remaining.erase(site);seeded=true;break;
    }
    if(!seeded && has_side) {current=side;remaining.erase(side);seeded=true;}
    if(seeded) {
        while(find_edge(graph,target,current)==graph.end()) {
            auto it=std::find_if(remaining.begin(),remaining.end(),[&](S s){return find_edge(graph,current,s)!=graph.end();});
            if(it==remaining.end())throw std::invalid_argument("disconnected native group");
            join(current,*it);
        }
        join(target,current);
    }
    while(!remaining.empty()) {
        auto it=std::find_if(remaining.begin(),remaining.end(),[&](S s){return find_edge(graph,target,s)!=graph.end();});
        if(it==remaining.end())throw std::invalid_argument("disconnected native group");
        join(target,*it);
    }
    return ops;
}

template<class S>
mcmc_edge_region_plan<S> mcmc_build_edge_region_plan(
    const std::vector<std::pair<S,S>>& original,const std::map<S,S>& representatives,
    const std::set<S>& targets,const std::map<S,S>& owners,const std::vector<S>& leaves,
    S target,bool has_side,S side,std::uint64_t count,std::uint64_t d,std::uint64_t chi,
    const std::vector<std::string>& edge_dimensions={}) {
    auto [edges,ports]=mcmc_edge_quotient(original,representatives);
    graph_string_shapes<S> shapes;std::map<S,int> ranks;
    for(auto item:representatives) {shapes[item.second];ranks[item.second]=static_cast<int>(owners.at(item.second));}
    for(auto& item:shapes)for(auto e:GetSurroundingBondIndex(item.first,edges)) {
        std::string dim;
        for(auto token:ports[e]) {
            const auto [a,b]=original[token];
            const auto term=edge_dimensions.empty()?std::string((!targets.count(a) && !targets.count(b))?"x":"D"):edge_dimensions.at(token);
            if(term!="1")dim+=term;
        }
        item.second.push_back(dim.empty()?"1":dim);
    }
    auto baseline=mcmc_edge_baseline(edges,leaves,target,has_side,side);
    const auto parsed=detail::parse_graph_symbolic_shapes(edges,shapes);
    auto choices=detail::build_graph_symbolic_region_plans_prepared(edges,parsed,leaves,ranks,target,count,baseline);
    const auto selected=select_graph_symbolic_region_plan(choices.first,choices.second,d,chi);
    auto ops=std::get<0>(selected<0?choices.second:choices.first.at(selected));
    // Independent common subtrees precede the repeated target spine. Re-evaluate
    // this valid topological order: its external fusion order is the one executed.
    std::set<S> dependent{target};std::vector<std::array<S,3>> common,repeated;
    for(auto st:ops) {
        if(dependent.count(st[0]) || dependent.count(st[1])) {
            repeated.push_back(st);dependent.erase(st[0]);dependent.erase(st[1]);dependent.insert(st[2]);
        }else common.push_back(st);
    }
    common.insert(common.end(),repeated.begin(),repeated.end());
    auto evaluated=detail::evaluate_graph_symbolic_region_prepared(edges,parsed,common,ranks,target,count);
    mcmc_edge_region_plan<S> result{std::move(edges),std::move(ports),std::move(common),std::get<3>(evaluated),target,{},0,{}, {}};
    mcmc_compile_edge_region(result);return result;
}

/** Symbolic shapes at the boundary-MPS checkpoint used by sampling. The
 * ordered graph is the result of mcmc_init_edges (including bond shifts).
 * Only an absorbed/truncated boundary has model chi on its internal bonds;
 * the outermost unabsorbed row still has original D bonds. This is metadata
 * only and deliberately does not predict the numerical truncated ranks.
 */
template<class S>
std::vector<std::string> mcmc_boundary_edge_model(const std::vector<std::pair<S,S>>& graph,
    const std::vector<std::vector<S>>& lines,std::size_t target) {
    if(target>=lines.size())throw std::invalid_argument("boundary model target line");
    std::set<S> lower,upper;
    if(target>1)lower.insert(lines[target-1].begin(),lines[target-1].end());
    if(target+2<lines.size())upper.insert(lines[target+1].begin(),lines[target+1].end());
    std::vector<std::string> dimensions;
    for(auto e:graph)dimensions.push_back((lower.count(e.first) && lower.count(e.second)) ||
        (upper.count(e.first) && upper.count(e.second))?"x":"D");
    return dimensions;
}

template<class S>
mcmc_edge_line_plan<S> mcmc_build_edge_line_plan(const std::vector<std::pair<S,S>>& edges,
    const std::vector<S>& line,const std::vector<std::vector<S>>& groups,const std::map<S,S>& owners,
    int direction,int width,std::uint64_t d=0,std::uint64_t chi=0,const std::map<S,int>& counts={},
    const std::vector<std::string>& edge_dimensions={}) {
    if((width!=1 && width!=2) || line.empty() || groups.size()!=line.size())throw std::invalid_argument("native line inputs");
    if(!edge_dimensions.empty() && edge_dimensions.size()!=edges.size())throw std::invalid_argument("native edge model size");
    for(const auto& dim:edge_dimensions)detail::parse_symbolic_axis(dim);
    std::map<S,S> identity;std::set<S> targets(line.begin(),line.end());
    for(std::size_t j=0;j<groups.size();++j) {
        if(groups[j].empty() || groups[j].front()!=line[j])throw std::invalid_argument("native group anchor");
        for(auto site:groups[j])if(!identity.emplace(site,site).second)throw std::invalid_argument("duplicate native group site");
    }
    for(auto e:edges)if(!identity.count(e.first) || !identity.count(e.second))throw std::invalid_argument("native edge outside groups");
    std::vector<std::size_t> order(line.size());std::iota(order.begin(),order.end(),0);if(direction)std::reverse(order.begin(),order.end());
    std::vector<mcmc_edge_region_plan<S>> suffix(line.size()+1);
    std::vector<mcmc_edge_window_plan<S>> windows;std::set<std::vector<int>> sets;
    auto add_set=[&](const auto& p){if(!std::get<3>(p).empty())sets.insert(std::get<3>(p));};
    auto collapse=[&](auto& reps,std::size_t first,std::size_t last,S root) {
        for(auto k=first;k<last;++k)for(auto site:groups[order[k]])reps[site]=root;
    };
    for(std::size_t end=line.size();end>static_cast<std::size_t>(width);) {
        --end;auto reps=identity;auto j=order[end];const bool side=end+1<line.size();
        const S side_site=side?line[order[end+1]]:line[j];
        if(side)collapse(reps,end+1,line.size(),side_site);
        auto leaves=groups[j];if(side)leaves.push_back(side_site);
        suffix[end]=mcmc_build_edge_region_plan(edges,reps,targets,owners,leaves,line[j],side,side_site,1,d,chi,edge_dimensions);add_set(suffix[end]);
    }
    for(std::size_t pos=0;pos+width<=line.size();++pos) {
        auto reps=identity;const auto a=order[pos],b=order[pos+width-1],end=pos+width;
        const bool left=pos>0,right=end<line.size();const S l=left?line[order[pos-1]]:line[a],r=right?line[order[end]]:line[b];
        if(left)collapse(reps,0,pos,l);if(right)collapse(reps,end,line.size(),r);
        auto leaves=groups[a];if(left)leaves.push_back(l);
        const auto p=counts.count(line[a])?counts.at(line[a]):1;
        auto ap=mcmc_build_edge_region_plan(edges,reps,targets,owners,leaves,line[a],left,l,p,d,chi,edge_dimensions);
        collapse(reps,0,pos+1,line[a]);
        mcmc_edge_region_plan<S> bp,jp;std::set<int> members(std::get<3>(ap).begin(),std::get<3>(ap).end());
        if(width==2) {
            auto other=groups[b];if(right)other.push_back(r);
            const auto q=counts.count(line[b])?counts.at(line[b]):1;
            bp=mcmc_build_edge_region_plan(edges,reps,targets,owners,other,line[b],right,r,q,d,chi,edge_dimensions);
            members.insert(std::get<3>(bp).begin(),std::get<3>(bp).end());collapse(reps,pos+1,line.size(),line[b]);
        }
        if(width==2 || right) {
            const auto partner=width==2?line[b]:r;
            auto [g,ports]=mcmc_edge_quotient(edges,reps);
            if(g.size()!=1 || find_edge(g,line[a],partner)==g.end())throw std::invalid_argument("window cut is not a single native edge");
            jp={g,ports,{{line[a],partner,line[a]}},{static_cast<int>(owners.at(line[a])),static_cast<int>(owners.at(partner))},line[a],{},0,{}, {}};
            mcmc_compile_edge_region(jp);
            members.insert(static_cast<int>(owners.at(partner)));
        }
        std::vector<int> participants(members.begin(),members.end());sets.insert(participants);
        windows.emplace_back(std::move(ap),std::move(bp),std::move(jp),std::move(participants),
            std::array<int,2>{counts.count(line[a])?counts.at(line[a]):0,width==2?(counts.count(line[b])?counts.at(line[b]):0):1});
    }
    return {std::move(order),std::move(suffix),std::move(windows),std::vector<std::vector<int>>(sets.begin(),sets.end())};
}
} // namespace tnbp
#endif
