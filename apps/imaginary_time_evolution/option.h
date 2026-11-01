#pragma once
#include <charconv>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>

namespace ite {
inline std::size_t count(const std::string& s) {
  std::size_t n=0;auto r=std::from_chars(s.data(),s.data()+s.size(),n);
  if(s.empty() || r.ec!=std::errc{} || r.ptr!=s.data()+s.size())
    throw std::invalid_argument("expected nonnegative integer: "+s);
  return n;
}
inline double real(const std::string& s) {
  std::size_t used=0;double n=std::stod(s,&used);
  if(used!=s.size() || !std::isfinite(n)) throw std::invalid_argument("expected finite real: "+s);
  return n;
}
struct Options {
  std::map<std::string,std::string> values;
  bool help=false;
  bool has(const std::string& key) const {return values.count(key);}
  std::string get(const std::string& key,const std::string& fallback="") const {
    auto it=values.find(key);return it==values.end() ? fallback : it->second;
  }
};
inline Options options(int argc,char** argv) {
  Options o;
  const std::set<std::string> valid={"inputedge","sparse_pauli","ownerfile","initial_state","initial_statefile","time_step","num_steps","schedule","trotter_order","max_bp_iterations","bp_tolerance","max_bond_dim","sv_min","truncation_error","measure_every","step_end","portable_save","portable_load","restart_mode"};
  for(int i=1;i<argc;++i) {
    const std::string flag=argv[i];
    if(flag=="--help" || flag=="-h") {o.help=true;continue;}
    if(flag.rfind("--",0)!=0 || !valid.count(flag.substr(2)))
      throw std::invalid_argument("unknown option: "+flag);
    if(i+1>=argc || std::string(argv[i+1]).empty() || std::string(argv[i+1]).rfind("--",0)==0)
      throw std::invalid_argument(flag+" requires a value");
    if(!o.values.emplace(flag.substr(2),argv[++i]).second) throw std::invalid_argument("duplicate option: "+flag);
  }
  if(o.has("initial_state") && o.has("initial_statefile")) throw std::invalid_argument("choose initial_state or initial_statefile");
  if(o.has("portable_load") && (o.has("initial_state") || o.has("initial_statefile"))) throw std::invalid_argument("portable_load conflicts with product state options");
  if(o.has("schedule") && (o.has("time_step") || o.has("num_steps"))) throw std::invalid_argument("choose schedule or time_step/num_steps");
  if(o.has("time_step")!=o.has("num_steps")) throw std::invalid_argument("time_step and num_steps are required together");
  if(o.get("restart_mode","auto")!="auto" && o.get("restart_mode")!="initial") throw std::invalid_argument("restart_mode must be auto or initial");
  if(o.has("restart_mode") && !o.has("portable_load")) throw std::invalid_argument("restart_mode requires portable_load");
  if(o.get("initial_state","zero")!="zero" && o.get("initial_state")!="plus") throw std::invalid_argument("initial_state must be zero or plus");
  return o;
}
} // namespace ite
