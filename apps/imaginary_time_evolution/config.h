#pragma once
#include "option.h"
#include "pauli/evolution_plan.h"
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>

namespace ite {
struct Stage { double time_step; std::size_t num_steps; };
struct Settings {
  std::size_t max_bp_iterations=50,max_bond_dim=100;
  double bp_tolerance=1e-8,sv_min=1e-8,truncation_error=1e-8;
  int trotter_order=2;
};
struct Config {tnbp::PauliHamiltonian h;std::vector<Stage> schedule;Settings settings;};
inline void require(bool ok,const std::string& message) {if(!ok)throw std::runtime_error("imaginary-time: "+message);}
inline std::size_t total_steps(const std::vector<Stage>& stages) {
  require(!stages.empty(),"empty schedule");std::size_t total=0;long double tau=0;
  for(auto s:stages) {
    require(std::isfinite(s.time_step) && s.time_step>0 && s.num_steps>0,"schedule requires positive finite time_step and positive num_steps");
    require(s.num_steps<=std::numeric_limits<std::size_t>::max()-total,"schedule step count overflow");total+=s.num_steps;
    tau+=static_cast<long double>(s.time_step)*s.num_steps;
    require(tau<=std::numeric_limits<double>::max(),"cumulative imaginary time overflow");
  }
  return total;
}
inline double elapsed_time(const std::vector<Stage>& stages,std::size_t next) {
  long double time=0;
  for(auto s:stages) {auto n=std::min(next,s.num_steps);time+=static_cast<long double>(s.time_step)*n;next-=n;}
  require(!next,"step exceeds schedule");return static_cast<double>(time);
}
inline double step_size(const std::vector<Stage>& stages,std::size_t step) {
  for(auto s:stages) {if(step<s.num_steps)return s.time_step;step-=s.num_steps;}
  throw std::invalid_argument("step outside schedule");
}
inline std::vector<Stage> read_schedule(const Options& o) {
  std::vector<Stage> stages;
  if(o.has("schedule")) {
    std::ifstream f(o.get("schedule"));require(bool(f),"cannot open schedule");std::string line;
    while(std::getline(f,line)) {
      line=line.substr(0,line.find('#'));std::istringstream row(line);std::string dt,n,extra;
      if(!(row>>dt))continue;
      require(bool(row>>n) && !(row>>extra),"schedule rows require time_step num_steps");
      stages.push_back({real(dt),count(n)});
    }
    require(!f.bad(),"schedule read failed");
  } else {
    require(o.has("time_step") && o.has("num_steps"),"specify time_step/num_steps or schedule");
    stages.push_back({real(o.get("time_step")),count(o.get("num_steps"))});
  }
  total_steps(stages);return stages;
}
inline void validate(const Settings& s) {
  require(s.max_bond_dim>0 && s.max_bond_dim<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"invalid max_bond_dim");
  require(std::isfinite(s.bp_tolerance) && s.bp_tolerance>=0 && std::isfinite(s.sv_min) && s.sv_min>=0 &&
    std::isfinite(s.truncation_error) && s.truncation_error>=0,"invalid BP/truncation tolerance");
  require(s.trotter_order==1 || s.trotter_order==2,"trotter_order must be 1 or 2");
}
inline Settings settings(const Options& o,Settings s={}) {
  if(o.has("max_bp_iterations"))s.max_bp_iterations=count(o.get("max_bp_iterations"));
  if(o.has("max_bond_dim"))s.max_bond_dim=count(o.get("max_bond_dim"));
  if(o.has("bp_tolerance"))s.bp_tolerance=real(o.get("bp_tolerance"));
  if(o.has("sv_min"))s.sv_min=real(o.get("sv_min"));
  if(o.has("truncation_error"))s.truncation_error=real(o.get("truncation_error"));
  if(o.has("trotter_order")) {auto order=count(o.get("trotter_order"));require(order==1 || order==2,"trotter_order must be 1 or 2");s.trotter_order=static_cast<int>(order);}
  validate(s);return s;
}
inline tnbp::PauliHamiltonian read_hamiltonian(const std::string& path) {
  // Existing SparsePauliOp text layout, with strict row validation.
  require(pauli::get_extension(path)=="txt","Hamiltonian extension must be .txt");
  std::ifstream f(path);require(bool(f),"cannot open Hamiltonian: "+path);
  std::vector<pauli::Term<std::complex<double>>> terms;std::string line;
  while(std::getline(f,line)) {
    line=line.substr(0,line.find('#'));auto slash=line.find("//");line=line.substr(0,slash);
    std::istringstream row(line);std::vector<std::string> tokens;std::string token;
    while(row>>token)tokens.push_back(token);
    if(tokens.empty())continue;
    require(tokens.size()==2 || tokens.size()==3,"Hamiltonian rows require real [imag] PauliString");
    terms.push_back({{real(tokens[0]),tokens.size()==3 ? real(tokens[1]):0},tokens.back()});
  }
  require(!f.bad(),"Hamiltonian read failed");return tnbp::NormalizePauliHamiltonian(terms);
}
inline std::string encode_config(const Config& c) {
  std::ostringstream out;out.imbue(std::locale::classic());out<<std::setprecision(17);
  out<<"N "<<c.h.num_sites<<" IDENTITY "<<c.h.identity_coefficient<<" TERMS "<<c.h.terms.size();
  for(const auto& t:c.h.terms)out<<' '<<t.coefficient<<' '<<t.label;
  out<<" SCHEDULE "<<c.schedule.size();for(auto s:c.schedule)out<<' '<<s.time_step<<' '<<s.num_steps;
  const auto& s=c.settings;
  out<<" SETTINGS "<<s.max_bp_iterations<<' '<<s.max_bond_dim<<' '<<s.bp_tolerance<<' '<<s.sv_min<<' '<<s.truncation_error<<' '<<s.trotter_order;
  return out.str();
}
inline Config decode_config(const std::string& text) {
  std::istringstream in(text);in.imbue(std::locale::classic());
  auto token=[&] {std::string s;require(bool(in>>s),"truncated configuration");return s;};
  auto tag=[&](const std::string& expected){require(token()==expected,"invalid configuration tag: "+expected);};
  tag("N");auto n=count(token());require(n>0 && n<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"invalid configuration site count");
  tag("IDENTITY");auto ci=real(token());tag("TERMS");auto k=count(token());require(k<=text.size(),"invalid configuration term count");
  std::vector<pauli::Term<std::complex<double>>> terms{{{ci,0},std::string(n,'I')}};
  for(std::size_t i=0;i<k;++i) {double c=real(token());auto label=token();require(label.size()==n,"configuration label length");terms.push_back({{c,0},label});}
  Config c;c.h=tnbp::NormalizePauliHamiltonian(terms);
  tag("SCHEDULE");k=count(token());require(k<=text.size(),"invalid schedule count");
  for(std::size_t i=0;i<k;++i) {double dt=real(token());auto steps=count(token());c.schedule.push_back({dt,steps});}
  total_steps(c.schedule);tag("SETTINGS");auto& s=c.settings;
  s.max_bp_iterations=count(token());s.max_bond_dim=count(token());s.bp_tolerance=real(token());s.sv_min=real(token());s.truncation_error=real(token());
  auto order=count(token());require(order==1 || order==2,"invalid Trotter order");s.trotter_order=static_cast<int>(order);validate(s);
  in>>std::ws;require(in.eof(),"trailing configuration data");return c;
}
} // namespace ite
