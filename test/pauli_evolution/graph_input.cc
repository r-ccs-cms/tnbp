#include "tnbp/framework/graphio.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
using Edges=std::vector<std::pair<int,int>>;
using Reader=Edges(*)(const std::string&,bool);
Reader legacy_reader=&tnbp::read_edges_from_file<int>;
int main(int argc,char**argv){
  if(argc!=2)return 2;
  const std::string root=argv[1],file=root+"/graph-input.txt";
  {std::ofstream out(file);out<<"2 1\n0 1\n1 2\n";}
  if(legacy_reader(file,true)!=Edges{{0,1},{1,2}})throw std::runtime_error("legacy canonical ordering");
  if(tnbp::read_edges_from_file<int>(file,false)!=Edges{{1,2},{0,1},{1,2}})throw std::runtime_error("input ordering");
  if(tnbp::read_edges_from_file<int>(file,true)!=Edges{{0,1},{1,2}})throw std::runtime_error("undirected deduplication");
  Edges normalized{{2,1},{0,1},{1,2}};
  tnbp::normalize_edges(normalized);
  if(normalized!=Edges{{1,2},{0,1},{1,2}})throw std::runtime_error("normalization preserves edge IDs");
  for(const auto& malformed:std::vector<std::string>{"oops 1\n","9999999999999999999999999 1\n"}) {
    {std::ofstream out(file);out<<malformed;}
    bool rejected=false;
    try{legacy_reader(file,false);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("malformed endpoint silently skipped");
  }
  {std::ofstream out(file);out<<"  # comment\n\n \t\n2 1 # edge\n";}
  if(legacy_reader(file,false)!=Edges{{1,2}})throw std::runtime_error("comments and blank rows");
  bool rejected=false;
  try{tnbp::read_edges_from_file<int>(root,false);}catch(const std::runtime_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("directory accepted as graph");
  std::cout<<"PASS graph input compatibility and read errors\n";
}
