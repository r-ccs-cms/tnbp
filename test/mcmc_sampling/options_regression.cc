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
  for(const char* flag:{"--initial_statefile","--portable_input","--inputname","--inputedge","--outputname","--seed",
      "--p_bond_dim","--v_bond_dim","--bmps_bond_dim","--bmps_bp_tolerance",
      "--bmps_max_bp_steps","--bmps_tg_err","--bmps_sv_min","--init_line",
      "--groupfile","--linefile","--num_samples","--num_batches","--sampling_method","--update_sites"}) {
    rejects({flag});rejects({flag,"--seed","1"});rejects({flag,""});
  }
  for(const char* flag:{"--p_bond_dim","--v_bond_dim","--bmps_bond_dim","--num_batches"})
    for(const char* value:{"0","-1","abc","1junk","2147483648","99999999999999999999999"}) rejects({flag,value});
  for(const char* flag:{"--num_samples","--bmps_max_bp_steps","--seed"})
    for(const char* value:{"-1","abc","1.5","1junk","99999999999999999999999"}) rejects({flag,value});
  for(const char* flag:{"--bmps_bp_tolerance","--bmps_tg_err","--bmps_sv_min"})
    for(const char* value:{"-1","nan","inf","abc","1e-8junk","1e9999"}) rejects({flag,value});
  for(const char* value:{"0,","0,,1",",0","x","2147483648"}) rejects({"--init_line",value});
  rejects({"--sampling_method","3"});
  rejects({"--update_sites","0"}); rejects({"--update_sites","3"});
  rejects({"--update_sites","2","--sampling_method","1"});
  if(parse({}).sampling_method!=2 || parse({}).update_sites!=1 ||
     parse({"--update_sites","2"}).update_sites!=2) throw std::runtime_error("full heat bath defaults");
  rejects({"--line_mode","other"});
  rejects({"--line_mode","alternating"});
  rejects({"--line_mode","alternating","--second_linefile","b","--update_sites","1"});
  rejects({"--line_mode","alternating","--second_linefile","b","--sampling_method","0"});
  for(const char* flag:{"--first_round_trips","--second_round_trips"}) {
    for(const char* value:{"0","-1","x","1.5","2147483648"}) rejects({flag,value});
    rejects({flag,"1"});
  }
  rejects({"--second_linefile","b"});rejects({"--second_groupfile","g"});
  const auto alternating=parse({"--line_mode","alternating","--second_linefile","b",
      "--second_groupfile","g","--update_sites","2","--first_round_trips","2","--second_round_trips","3"});
  if(alternating.first_round_trips!=2 || alternating.second_round_trips!=3 ||
     alternating.second_groupfile!="g" || alternating.line_mode!="alternating")
    throw std::runtime_error("alternating options changed");
  if(parse({}).line_mode!="single") throw std::runtime_error("single layout default changed");
  if(parse({"--initial_statefile","states.txt"}).initial_statefile!="states.txt")
    throw std::runtime_error("initial state option changed");
  rejects({"--unknown","1"});
  const auto result=parse({"--seed","4294967295","--num_samples","0","--bmps_max_bp_steps","0",
      "--num_batches","2","--init_line","-2,0,4","--bmps_sv_min","0","--sampling_method","1"});
  if(result.seed!=4294967295u || result.num_samples!=0 || result.bmps_max_bp_steps!=0 ||
     result.num_batches!=2 || result.init_line!=std::vector<int>{-2,0,4} ||
     result.bmps_sv_min!=0 || result.sampling_method!=1) throw std::runtime_error("valid options changed");
  if(!parse({"--num_samples","--help"}).help || !parse({"-h"}).help)
    throw std::runtime_error("help precedence changed");
  std::cout << "options regression PASS\n";
}
