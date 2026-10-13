using Real=double;
#include "../../apps/mcmc_sampling/option.h"
#include <functional>

Option parse(std::vector<std::string> args) {
  args.insert(args.begin(),"sampler");
  std::vector<char*> argv;
  for(auto& arg:args) argv.push_back(arg.data());
  return generate_options(static_cast<int>(argv.size()),argv.data());
}
void rejects(const std::vector<std::string>& args) {
  try {parse(args);}
  catch(const std::invalid_argument&) {return;}
  throw std::runtime_error("invalid option was accepted: "+args.front());
}
int main() {
  for(const char* flag:{"--inputname","--inputedge","--outputname","--seed",
      "--p_bond_dim","--v_bond_dim","--bmps_bond_dim","--bmps_bp_tolerance",
      "--bmps_max_bp_steps","--bmps_tg_err","--bmps_sv_min","--init_line",
      "--groupfile","--linefile","--num_samples","--num_batches","--sampling_method"}) {
    rejects({flag});rejects({flag,"--seed","1"});rejects({flag,""});
  }
  for(const char* flag:{"--p_bond_dim","--v_bond_dim","--bmps_bond_dim","--num_batches"})
    for(const char* value:{"0","-1","abc","1junk","2147483648","99999999999999999999999"}) rejects({flag,value});
  for(const char* flag:{"--num_samples","--bmps_max_bp_steps","--seed"})
    for(const char* value:{"-1","abc","1.5","1junk","99999999999999999999999"}) rejects({flag,value});
  for(const char* flag:{"--bmps_bp_tolerance","--bmps_tg_err","--bmps_sv_min"})
    for(const char* value:{"-1","nan","inf","abc","1e-8junk","1e9999"}) rejects({flag,value});
  for(const char* value:{"0,","0,,1",",0","x","2147483648"}) rejects({"--init_line",value});
  rejects({"--sampling_method","2"});rejects({"--unknown","1"});
  const auto result=parse({"--seed","4294967295","--num_samples","0","--bmps_max_bp_steps","0",
      "--num_batches","2","--init_line","-2,0,4","--bmps_sv_min","0","--sampling_method","1"});
  if(result.seed!=4294967295u || result.num_samples!=0 || result.bmps_max_bp_steps!=0 ||
     result.num_batches!=2 || result.init_line!=std::vector<int>{-2,0,4} ||
     result.bmps_sv_min!=0 || result.sampling_method!=1) throw std::runtime_error("valid options changed");
  if(!parse({"--num_samples","--help"}).help || !parse({"-h"}).help)
    throw std::runtime_error("help precedence changed");
  std::cout << "options regression PASS\n";
}
