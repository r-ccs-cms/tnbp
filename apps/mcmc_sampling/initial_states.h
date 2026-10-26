#pragma once
#include <charconv>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Return [global site position][batch] in the supplied site order. File order
// is arbitrary; all sites (including fixed ancillas) must be specified once.
inline std::vector<std::vector<int>> read_initial_states(
    const std::string& filename,const std::vector<int>& sites,std::size_t batches) {
  std::ifstream input(filename);
  if(!input) throw std::runtime_error("cannot open initial state file: "+filename);
  std::map<int,std::size_t> positions;
  for(std::size_t i=0;i<sites.size();++i) positions.emplace(sites[i],i);
  std::vector<std::vector<int>> states(sites.size());
  std::vector<bool> seen(sites.size(),false);
  std::string line;
  std::size_t line_number=0;
  while(std::getline(input,line)) {
    ++line_number;
    line.erase(line.find('#')==std::string::npos ? line.size() : line.find('#'));
    std::istringstream row(line);
    std::vector<int> values;
    std::string token;
    while(row>>token) {
      int value;
      auto result=std::from_chars(token.data(),token.data()+token.size(),value);
      if(result.ec!=std::errc{} || result.ptr!=token.data()+token.size())
        throw std::invalid_argument("initial state file: invalid integer on line "+std::to_string(line_number));
      values.push_back(value);
    }
    if(values.empty()) continue;
    if(values.size()!=batches+1)
      throw std::invalid_argument("initial state file: expected site ID and one state per batch on line "+std::to_string(line_number));
    auto found=positions.find(values[0]);
    if(found==positions.end()) throw std::invalid_argument("initial state file: unknown site "+std::to_string(values[0]));
    const auto i=found->second;
    if(seen[i]) throw std::invalid_argument("initial state file: duplicate site "+std::to_string(values[0]));
    seen[i]=true;
    states[i].assign(values.begin()+1,values.end());
    for(int state:states[i]) if(state<0)
      throw std::invalid_argument("initial state file: negative physical state");
  }
  if(input.bad()) throw std::runtime_error("cannot read initial state file: "+filename);
  for(std::size_t i=0;i<sites.size();++i) if(!seen[i])
    throw std::invalid_argument("initial state file: missing site "+std::to_string(sites[i]));
  return states;
}

// Only owners need tensor dimensions; no global dimension exchange is needed.
inline std::vector<std::vector<int>> local_initial_states(
    const std::vector<std::vector<int>>& states,const std::vector<int>& sites,
    const std::vector<int>& local_sites,const std::vector<int>& dimensions,std::size_t batches) {
  std::map<int,std::size_t> positions;
  for(std::size_t i=0;i<sites.size();++i) positions.emplace(sites[i],i);
  std::vector<std::vector<int>> local(batches,std::vector<int>(local_sites.size()));
  for(std::size_t i=0;i<local_sites.size();++i) {
    const auto& row=states.at(positions.at(local_sites[i]));
    for(std::size_t batch=0;batch<batches;++batch) {
      const auto state=row.at(batch);
      if(state<0 || state>=dimensions.at(i))
        throw std::invalid_argument("initial state file: physical state out of range at site "+std::to_string(local_sites[i]));
      local[batch][i]=state;
    }
  }
  return local;
}
