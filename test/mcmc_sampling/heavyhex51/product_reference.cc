// Independent scalar Metropolis / heat-bath simulation for the D=1 probe, no TNBP/TCAPI.
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <cmath>
int main(int argc,char**argv) {
 if(argc!=2) return 2;
 std::ifstream f(argv[1]);std::string row;std::vector<std::vector<int>> lines;std::vector<std::string> actual;
 int ranks=1,method=0;
 bool uniform=std::string(argv[1]).find("uniform-")!=std::string::npos;
 while(std::getline(f,row)) {
  if(row.rfind("mode=",0)==0) {
   ranks=std::stoi(row.substr(row.find("ranks=")+6));
   if(row.find("method=")!=std::string::npos) method=std::stoi(row.substr(row.find("method=")+7));
   if(row.find("uniform=")!=std::string::npos) uniform=std::stoi(row.substr(row.find("uniform=")+8));
  }
  if(row.rfind("line ",0)==0) {std::istringstream in(row.substr(row.find(':')+1));std::vector<int> l;int v;while(in>>v)l.push_back(v);lines.push_back(l);}
  if(row.rfind("sample ",0)==0) actual.push_back(row.substr(row.find(':')+2));
 }
 std::vector<int> sites;for(auto l:lines)sites.insert(sites.end(),l.begin(),l.end());std::sort(sites.begin(),sites.end());
 std::map<int,int> state;std::map<int,double> probability;std::map<int,int> owners;
 for(std::size_t j=0;j<sites.size();++j) {
  double p=uniform ? .5 : .25+.08*(j%5),beta=uniform ? 0 : .37;
  double a=std::sin(beta/2)*std::sqrt(1-p)+std::cos(beta/2)*std::sqrt(p);
  probability[sites[j]]=a*a;state[sites[j]]=0;owners[sites[j]]=j%ranks;
 }
 std::vector<std::mt19937> generators;for(int r=0;r<ranks;++r) generators.emplace_back(7129+r);std::uniform_int_distribution<int> proposal(0,0);std::uniform_real_distribution<double> accept(0,1);
 auto step=[&](int v){auto& rng=generators[owners.at(v)];proposal(rng);const double p=probability.at(v),ratio=state[v] ? (1-p)/p : p/(1-p);if(accept(rng)<(method==1 ? (state[v] ? 1-p : p) : ratio))state[v]=1-state[v];};
 for(std::size_t n=0;n<16+actual.size();++n) {
  for(const auto& l:lines)for(auto v:l)step(v);
  for(auto l=lines.rbegin();l!=lines.rend();++l)for(auto v=l->rbegin();v!=l->rend();++v)step(*v);
  if(n>=16) {
   std::string expected;for(auto v:sites)expected+=char('0'+state[v]);
   if(expected!=actual[n-16]) {std::cerr<<"MISMATCH sample "<<n-16<<'\n';return 1;}
  }
 }
 if(actual.empty()) return 2;
 std::cout<<"PASS independent scalar sampler: "<<actual.size()<<" complete bitstrings match\n";
}
