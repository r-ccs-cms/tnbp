#ifndef TNBP_PTNS_FUNCTION_MCMC_PLAN_H
#define TNBP_PTNS_FUNCTION_MCMC_PLAN_H
#include "tnbp/ptns/function/bmps.h"
#include <set>
#include <stdexcept>
#include <string>

namespace tnbp {
// groups[line][position]: target site first, followed by neighboring-line sites.
// Member order is the explicit absorption order; directions are specified separately.
template<class IntT>
using MCMCGrouping = std::vector<std::vector<std::vector<IntT>>>;

// Shared structural contract for automatic and manual sampling plans.
template<class IntT>
std::map<IntT,std::size_t> validate_mcmc_lines(
    const std::vector<std::pair<IntT,IntT>>& edges,
    const std::vector<std::vector<IntT>>& lines) {
    if(lines.size()<2)
        throw std::invalid_argument("sampling requires at least two nonempty lines");
    std::map<IntT,std::size_t> layer;
    for(std::size_t l=0;l<lines.size();++l) {
        if(lines[l].empty()) throw std::invalid_argument("empty line");
        for(auto s:lines[l]) if(!layer.emplace(s,l).second)
            throw std::invalid_argument("duplicate site in lines");
    }
    const auto sites=GetSiteIndexFromBond(edges);
    if(sites.size()!=layer.size()) throw std::invalid_argument("lines must cover exactly all graph sites");
    for(auto s:sites) if(!layer.count(s)) throw std::invalid_argument("graph site missing from lines");
    for(auto [a,b]:edges) {
        auto x=layer.at(a),y=layer.at(b);
        if(a==b || (x>y ? x-y : y-x)>1)
            throw std::invalid_argument("sampling supports edges within or between adjacent lines only");
    }
    return layer;
}

template<class IntT>
void validate_mcmc_grouping(const std::vector<std::pair<IntT,IntT>>& edges,
    const std::vector<std::vector<IntT>>& lines,const MCMCGrouping<IntT>& groups) {
    const auto layer=validate_mcmc_lines(edges,lines);
    if(groups.size()!=lines.size())
        throw std::invalid_argument("manual grouping requires one group list per line");
    for(std::size_t l=0;l<lines.size();++l) {
        if(groups[l].size()!=lines[l].size()) throw std::invalid_argument("missing target group");
        std::set<IntT> expected(lines[l].begin(),lines[l].end()),seen;
        if(l) expected.insert(lines[l-1].begin(),lines[l-1].end());
        if(l+1<lines.size()) expected.insert(lines[l+1].begin(),lines[l+1].end());
        for(std::size_t j=0;j<lines[l].size();++j) {
            const auto& g=groups[l][j];
            if(g.empty() || g.front()!=lines[l][j]) throw std::invalid_argument("wrong target group anchor");
            for(std::size_t k=0;k<g.size();++k) {
                if(!expected.count(g[k]) || !seen.insert(g[k]).second)
                    throw std::invalid_argument("unknown, nonadjacent or duplicate group member");
                if(k && layer.at(g[k])==l) throw std::invalid_argument("group contains another target-line site");
            }
        }
        if(seen!=expected) throw std::invalid_argument("group partition omits neighboring sites");
    }
    // Validate both directions on the evolving graph after absorption AND shifts.
    auto lower=edges,upper=edges;
    for(std::size_t l=1;l<lines.size();++l)
        bmps_increment_line_for_edges(lines[l-1],lines[l],lower,&groups[l]);
    for(std::size_t l=lines.size()-1;l>0;--l)
        bmps_increment_line_for_edges(lines[l],lines[l-1],upper,&groups[l-1]);
    // Each environment group and each growing left/right environment must be connected.
    for(std::size_t l=0;l<lines.size();++l) {
        auto local=edges;
        for(std::size_t k=1;k<l;++k) bmps_increment_line_for_edges(lines[k-1],lines[k],local,&groups[k]);
        for(std::size_t k=lines.size()-1;k>l+1;--k) bmps_increment_line_for_edges(lines[k],lines[k-1],local,&groups[k-1]);
        auto connected=[&](const std::set<IntT>& members) {
            std::set<IntT> reached{*members.begin()};
            bool changed=true;
            while(changed) { changed=false;
                for(auto [a,b]:local) if(members.count(a)&&members.count(b)) {
                    if(reached.count(a)) changed=reached.insert(b).second||changed;
                    if(reached.count(b)) changed=reached.insert(a).second||changed;
                }
            }
            if(reached!=members) throw std::invalid_argument("disconnected environment group or sweep prefix/suffix");
        };
        std::set<IntT> left,right;
        for(std::size_t j=0;j<groups[l].size();++j) {
            connected(std::set<IntT>(groups[l][j].begin(),groups[l][j].end()));
            left.insert(groups[l][j].begin(),groups[l][j].end());connected(left);
            const auto& g=groups[l][groups[l].size()-1-j];
            right.insert(g.begin(),g.end());connected(right);
        }
    }
}
// Build once from the original graph, independently for each consecutive
// three-line window. Lower members precede upper members in each anchored group.
// No tensor contraction or mutation of the original edge order is involved.
template<class IntT>
MCMCGrouping<IntT> make_mcmc_grouping(
    const std::vector<std::pair<IntT,IntT>>& edges,
    const std::vector<std::vector<IntT>>& lines) {
    validate_mcmc_lines(edges,lines);
    MCMCGrouping<IntT> groups(lines.size());
    for(std::size_t l=0;l<lines.size();++l) {
        for(auto site:lines[l]) groups[l].push_back({site});
        for(int direction:{-1,1}) {
            if((direction<0 && l==0) || (direction>0 && l+1==lines.size())) continue;
            const auto adjacent=direction<0 ? l-1 : l+1;
            const auto side=bmps_neighbor_groups(edges,lines[l],lines[adjacent]);
            for(std::size_t j=0;j<side.size();++j)
                groups[l][j].insert(groups[l][j].end(),side[j].begin()+1,side[j].end());
        }
    }
    // A nearest-neighbor greedy partition is not a solver for arbitrary line
    // layouts. Reject incompatible absorption/sweep connectivity before mutation.
    try { validate_mcmc_grouping(edges,lines,groups); }
    catch(const std::invalid_argument& e) {
        throw std::invalid_argument(std::string("automatic grouping is incompatible with this line layout: ")+
                                    e.what()+"; specify manual groups or change the lines");
    }
    return groups;
}
} // namespace tnbp
#endif
