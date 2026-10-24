#pragma once
#include "tnbp/framework/mpiutility.h"
#include <algorithm>
#include <charconv>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// CLI policy: ownership files and logging belong to the applications.
namespace app_ownership {
inline std::map<int,int> select(const std::vector<int>& sites,
                               const std::string& filename, MPI_Comm comm,
                               bool portable=false) {
  int rank,size; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
  std::map<int,int> owners;
  if(filename.empty()) {
    if(portable && sites.size()<static_cast<std::size_t>(size))
      throw std::invalid_argument("empty ranks are unsupported for portable input");
    return tnbp::ContiguousSiteOwners(sites,comm);
  }
  std::string error;
  std::vector<int> ranks(sites.size());
  if(rank==0) try {
    std::ifstream file(filename);
    if(!file) throw std::runtime_error("cannot open ownerfile: "+filename);
    std::string line; unsigned row=0;
    auto integer=[](const std::string& token) {
      int value=0; auto parsed=std::from_chars(token.data(),token.data()+token.size(),value);
      if(parsed.ec!=std::errc{} || parsed.ptr!=token.data()+token.size())
        throw std::runtime_error("ownerfile requires integer site_id mpi_rank");
      return value;
    };
    while(std::getline(file,line)) {
      ++row; line=line.substr(0,line.find('#')); std::istringstream in(line);
      std::string site_text,rank_text,extra;
      if(!(in>>site_text)) continue;
      if(!(in>>rank_text) || (in>>extra))
        throw std::runtime_error("ownerfile row "+std::to_string(row)+" requires site_id mpi_rank");
      int site=integer(site_text),owner=integer(rank_text);
      if(!std::binary_search(sites.begin(),sites.end(),site))
        throw std::runtime_error("unknown site in ownerfile: "+site_text);
      if(owner<0 || owner>=size) throw std::runtime_error("ownerfile rank out of range");
      if(!owners.emplace(site,owner).second) throw std::runtime_error("duplicate site in ownerfile");
    }
    if(file.bad()) throw std::runtime_error("failed reading ownerfile");
    if(owners.size()!=sites.size()) throw std::runtime_error("missing site in ownerfile");
    std::vector<int> count(size);
    for(std::size_t i=0;i<sites.size();++i) { ranks[i]=owners.at(sites[i]); ++count[ranks[i]]; }
    if(portable && std::find(count.begin(),count.end(),0)!=count.end())
      throw std::runtime_error("portable ownership does not support empty ranks");
  } catch(const std::exception& e) { error=e.what(); }
  tnbp::MpiBcast(error,0,comm);
  if(!error.empty()) throw std::runtime_error(error);
  MPI_Bcast(ranks.data(),static_cast<int>(ranks.size()),MPI_INT,0,comm);
  owners.clear();
  for(std::size_t i=0;i<sites.size();++i) owners.emplace(sites[i],ranks[i]);
  return owners;
}
inline void print(const std::map<int,int>& owners, const std::string& source, MPI_Comm comm) {
  int rank,size; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
  if(rank!=0) return;
  std::cout<<"# site ownership: "<<source<<'\n';
  for(int r=0;r<size;++r) {
    std::cout<<"# MPI rank "<<r<<":";
    for(const auto& [site,owner]:owners) if(owner==r) std::cout<<' '<<site;
    std::cout<<'\n';
  }
}
inline void check_legacy(const std::map<int,int>& stored,const std::map<int,int>& requested,
                         MPI_Comm comm) {
  int mismatch=stored!=requested, any=0;
  MPI_Allreduce(&mismatch,&any,1,MPI_INT,MPI_MAX,comm);
  if(any) throw std::runtime_error("ownerfile differs from legacy checkpoint ownership; use portable input to repartition");
}
} // namespace app_ownership
