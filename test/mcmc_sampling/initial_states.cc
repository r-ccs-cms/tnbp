#include "../../apps/mcmc_sampling/initial_states.h"
#include <iostream>
int main(int argc,char** argv) {
  if(argc!=2) return 1;
  const std::string file=argv[1];
  std::ofstream(file)<<"# site batch0 batch1\n20 5 6 # comment\n\n7 3 12\n";
  const auto states=read_initial_states(file,{7,20},2);
  if(local_initial_states(states,{7,20},{20,7},{16,16},2)!=
     std::vector<std::vector<int>>{{5,3},{6,12}})
    throw std::runtime_error("local order/batch mapping mismatch");
  if(local_initial_states(states,{7,20},{},{},2)!=std::vector<std::vector<int>>(2))
    throw std::runtime_error("empty owner mismatch");
  // A dimension-one ancilla accepts only zero; physical dimensions are per-site.
  const std::vector<std::vector<int>> mixed{{0,0},{3,12}};
  if(local_initial_states(mixed,{7,20},{7,20},{1,16},2)!=
     std::vector<std::vector<int>>{{0,3},{0,12}}) throw std::runtime_error("mixed dimensions");
  bool rejected=false;
  try {local_initial_states(mixed,{7,20},{20},{4},2);}
  catch(const std::invalid_argument&) {rejected=true;}
  if(!rejected) throw std::runtime_error("out of range accepted");
  std::cout<<"initial state local mapping PASS\n";
}
