#ifndef TNBP_PTNS_FUNCTION_MCMC_GROUP_IO_H
#define TNBP_PTNS_FUNCTION_MCMC_GROUP_IO_H
#include "tnbp/ptns/function/mcmc_plan.h"
#include <fstream>
#include <sstream>

namespace tnbp {
// Each non-comment row: zero-based line index, target site, ordered members.
template<class IntT>
MCMCGrouping<IntT> read_mcmc_grouping(std::istream& input,
                                     const std::vector<std::vector<IntT>>& lines) {
    MCMCGrouping<IntT> result(lines.size());
    for(std::size_t l=0;l<lines.size();++l) result[l].resize(lines[l].size());
    std::string row;std::size_t number=0;
    while(std::getline(input,row)) {
        ++number;row=row.substr(0,row.find('#'));
        std::istringstream stream(row);
        stream>>std::ws;if(stream.eof()) continue;
        long long line;IntT target;
        if(!(stream>>line>>target) || line<0 || static_cast<std::size_t>(line)>=lines.size())
            throw std::invalid_argument("invalid grouping row "+std::to_string(number));
        auto it=std::find(lines[line].begin(),lines[line].end(),target);
        if(it==lines[line].end()) throw std::invalid_argument("group target not in indicated line");
        auto& group=result[line][std::distance(lines[line].begin(),it)];
        if(!group.empty()) throw std::invalid_argument("duplicate target group row");
        group.push_back(target);
        while(true) {
            stream>>std::ws;if(stream.eof()) break;
            IntT member;
            if(!(stream>>member)) throw std::invalid_argument("invalid group member on row "+std::to_string(number));
            group.push_back(member);
        }
    }
    for(const auto& line:result) for(const auto& group:line)
        if(group.empty()) throw std::invalid_argument("missing target group row");
    return result;
}
}
#endif
