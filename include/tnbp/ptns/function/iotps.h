/// This file is a part of r-ccs-cms/tnbp
/**
   @file tnbp/ptns/function/bp.h
   @brief Belief propagation methods
*/
#ifndef TNBP_PTNS_FUNCTION_IOTPS_H
#define TNBP_PTNS_FUNCTION_IOTPS_H

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <iomanip>
#include <type_traits>

namespace tnbp {

  template <typename IntT>
  std::string rank_to_string(IntT rank, int width = 6) {
    static_assert(std::is_integral<IntT>::value, "IntT must be an integral type.");
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << rank;
    return oss.str();
  }
  
  template <typename IntT>
  std::string make_rank_filename(const std::string & prefix, IntT rank,
				 const std::string & ext = ".dat",
				 int width = 6) {
    return prefix + rank_to_string(rank,width) + ext;
  }
			       
  
  template <typename TenT>
  void SaveTPS(tci::context_handle_t<TenT> & ctx,
	       std::ostream & out,
	       const std::vector<TenT> & V,
	       const std::vector<int> & SiteIdx,
	       const std::map<int,int> & Site_To_MpiRank,
	       const std::vector<TenT> & E,
	     const std::vector<int> & EdgeIdx) {
    size_t size_V = V.size();
    size_t size_M = Site_To_MpiRank.size();
    size_t size_E = EdgeIdx.size();
    out.write(reinterpret_cast<char *>(&size_V),sizeof(size_t));
    out.write(reinterpret_cast<char *>(&size_M),sizeof(size_t));
    out.write(reinterpret_cast<char *>(&size_E),sizeof(size_t));
    for(const auto & Vi : V) {
      tci::save(ctx,Vi,out);
    }
    for(const auto & Idx : SiteIdx) {
      out.write(reinterpret_cast<const char *>(&Idx),sizeof(int));
    }
    for(const auto & [key,value] : Site_To_MpiRank) {
      out.write(reinterpret_cast<const char*>(&key), sizeof(int));
      out.write(reinterpret_cast<const char*>(&value), sizeof(int));
    }
    for(const auto & Em : E) {
      tci::save(ctx,Em,out);
    }
    for(const auto & Idx : EdgeIdx) {
      out.write(reinterpret_cast<const char *>(&Idx),sizeof(int));
    }
  }
  
  template <typename TenT>
  void LoadTPS(tci::context_handle_t<TenT> & ctx,
	       std::istream & in,
	       std::vector<TenT> & V,
	       std::vector<int> & SiteIdx,
	       std::map<int,int> & Site_To_MpiRank,
	       std::vector<TenT> & E,
	       std::vector<int> & EdgeIdx) {
    size_t size_V;
    size_t size_M;
    size_t size_E;
    in.read(reinterpret_cast<char *>(&size_V),sizeof(size_t));
    in.read(reinterpret_cast<char *>(&size_M),sizeof(size_t));
    in.read(reinterpret_cast<char *>(&size_E),sizeof(size_t));
    V.resize(size_V);
    SiteIdx.resize(size_V);
    E.resize(2*size_E);
    EdgeIdx.resize(size_E);
    for(size_t i=0; i < size_V; i++) {
      V[i] = tci::load<TenT>(ctx,in);
    }
    for(size_t i=0; i < size_V; i++) {
      in.read(reinterpret_cast<char *>(&SiteIdx[i]),sizeof(int));
    }
    for(size_t k=0; k < size_M; k++) {
      int key, value;
      in.read(reinterpret_cast<char *>(&key),sizeof(int));
      in.read(reinterpret_cast<char *>(&value),sizeof(int));
      Site_To_MpiRank[key] = value;
    }
    for(size_t m=0; m < 2*size_E; m++) {
      E[m] = tci::load<TenT>(ctx,in);
    }
    for(size_t m=0; m < size_E; m++) {
      in.read(reinterpret_cast<char *>(&EdgeIdx[m]),sizeof(int));
    }
  }

}

#endif
