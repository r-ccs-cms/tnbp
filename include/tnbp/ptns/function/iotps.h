/// This file is a part of r-ccs-cms/tnbp
/**
   @file tnbp/ptns/function/bp.h
   @brief Belief propagation methods
*/
#ifndef TNBP_PTNS_FUNCTION_IOTPS_H
#define TNBP_PTNS_FUNCTION_IOTPS_H

#include <istream>
#include <cstdint>
#include "tnbp/ptns/function/message_gauge.h"

#include "tnbp/framework/typedef.h"
#include <cstddef>
#include <map>

#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <iomanip>
#include <type_traits>

namespace tnbp {

  struct BPCheckpointMetadata {
    std::string source_algorithm = "unspecified";
    std::string source_checkpoint;
    std::string source_fingerprint = "unavailable";
  };

  // Non-cryptographic byte fingerprint for provenance. Runners that require
  // integrity/security guarantees should additionally retain a SHA-256 manifest.
  inline std::string CheckpointFingerprint(std::istream& in) {
    const auto state=in.rdstate();
    const auto start=in.tellg();
    if (start==std::istream::pos_type(-1)) { in.clear(state); return "unavailable"; }
    std::uint64_t hash=14695981039346656037ULL;
    char buffer[8192];
    while (in.read(buffer,sizeof(buffer)) || in.gcount())
      for(std::streamsize i=0;i<in.gcount();++i) {
        hash^=static_cast<unsigned char>(buffer[i]); hash*=1099511628211ULL;
      }
    in.clear(); in.seekg(start);
    if(!in) throw std::runtime_error("cannot restore checkpoint stream after fingerprint");
    in.clear(state);
    std::ostringstream out; out << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
  }

  // A trailing, versioned record preserves the existing tensor payload layout.
  // Legacy readers stop before this record; legacy files remain loadable.
  inline void WriteBPMetadata(std::ostream& out, const BPCheckpointMetadata& metadata) {
    out << "TNBP_BP_METADATA_V1\n" << std::quoted(BPAlgorithm) << '\n'
        << std::quoted(BPGauge) << '\n' << std::quoted(BPResidual) << '\n'
        << std::quoted(metadata.source_algorithm) << '\n'
        << std::quoted(metadata.source_checkpoint) << '\n'
        << std::quoted(metadata.source_fingerprint) << '\n';
    if (!out) throw std::runtime_error("failed to write BP checkpoint metadata");
  }


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
  void SaveTPS(tcapi::context_handle_t<TenT> & ctx,
	       std::ostream & out,
	       const std::vector<TenT> & V,
	       const std::vector<int> & SiteIdx,
	       const std::map<int,int> & Site_To_MpiRank,
	       const std::vector<TenT> & E,
	     const std::vector<int> & EdgeIdx,
               const BPCheckpointMetadata& metadata = {}) {
    std::size_t size_V = V.size();
    std::size_t size_M = Site_To_MpiRank.size();
    std::size_t size_E = EdgeIdx.size();
    out.write(reinterpret_cast<char *>(&size_V),sizeof(std::size_t));
    out.write(reinterpret_cast<char *>(&size_M),sizeof(std::size_t));
    out.write(reinterpret_cast<char *>(&size_E),sizeof(std::size_t));
    for(const auto & Vi : V) {
      tcapi::save(ctx,Vi,out);
    }
    for(const auto & Idx : SiteIdx) {
      out.write(reinterpret_cast<const char *>(&Idx),sizeof(int));
    }
    for(const auto & [key,value] : Site_To_MpiRank) {
      out.write(reinterpret_cast<const char*>(&key), sizeof(int));
      out.write(reinterpret_cast<const char*>(&value), sizeof(int));
    }
    for(std::size_t m=0;m<E.size();++m) {
      auto canonical=tcapi::copy(ctx,E[m]);
      FixMessageTracePhase(ctx,canonical,"save message="+std::to_string(m));
      tcapi::save(ctx,canonical,out);
    }
    for(const auto & Idx : EdgeIdx) {
      out.write(reinterpret_cast<const char *>(&Idx),sizeof(int));
    }
    WriteBPMetadata(out,metadata);
  }
  
  template <typename TenT>
  void LoadTPS(tcapi::context_handle_t<TenT> & ctx,
	       std::istream & in,
	       std::vector<TenT> & V,
	       std::vector<int> & SiteIdx,
	       std::map<int,int> & Site_To_MpiRank,
	       std::vector<TenT> & E,
	       std::vector<int> & EdgeIdx,
               BPCheckpointMetadata* metadata = nullptr) {
    const auto fingerprint=CheckpointFingerprint(in);
    std::size_t size_V=0;
    std::size_t size_M=0;
    std::size_t size_E=0;
    in.read(reinterpret_cast<char *>(&size_V),sizeof(std::size_t));
    in.read(reinterpret_cast<char *>(&size_M),sizeof(std::size_t));
    in.read(reinterpret_cast<char *>(&size_E),sizeof(std::size_t));
    if (!in) throw std::runtime_error("truncated TPS checkpoint header");
    Site_To_MpiRank.clear();
    V.resize(size_V);
    SiteIdx.resize(size_V);
    E.resize(2*size_E);
    EdgeIdx.resize(size_E);
    for(std::size_t i=0; i < size_V; i++) {
      V[i] = tcapi::load<TenT>(ctx,in);
    }
    for(std::size_t i=0; i < size_V; i++) {
      in.read(reinterpret_cast<char *>(&SiteIdx[i]),sizeof(int));
    }
    for(std::size_t k=0; k < size_M; k++) {
      int key, value;
      in.read(reinterpret_cast<char *>(&key),sizeof(int));
      in.read(reinterpret_cast<char *>(&value),sizeof(int));
      Site_To_MpiRank[key] = value;
    }
    for(std::size_t m=0; m < 2*size_E; m++) {
      E[m] = tcapi::load<TenT>(ctx,in);
    }
    for(std::size_t m=0; m < size_E; m++) {
      in.read(reinterpret_cast<char *>(&EdgeIdx[m]),sizeof(int));
    }
    if (!in) throw std::runtime_error("truncated TPS checkpoint");
    std::string algorithm="legacy-unspecified";
    if (in.peek()!=std::char_traits<char>::eof()) {
      std::string marker,gauge,residual,parent,source,parent_fingerprint;
      std::getline(in,marker);
      if (marker!="TNBP_BP_METADATA_V1" ||
          !(in >> std::quoted(algorithm) >> std::quoted(gauge) >> std::quoted(residual)
               >> std::quoted(parent) >> std::quoted(source) >> std::quoted(parent_fingerprint)))
        throw std::runtime_error("invalid BP checkpoint metadata");
      in.get(); // final newline
      if (algorithm!=BPAlgorithm || gauge!=BPGauge || residual!=BPResidual)
        throw std::runtime_error("unsupported BP checkpoint algorithm or convention");
    } else {
      in.clear(in.rdstate() & ~std::ios::eofbit);
    }
    NormalizeMessages(ctx,E);
    if (metadata) {
      metadata->source_algorithm=algorithm;
      metadata->source_fingerprint=fingerprint;
    }
    std::clog << "BP_ALGORITHM " << BPAlgorithm << " gauge=" << BPGauge
              << " residual=" << BPResidual << " source_algorithm=" << algorithm
              << " source_fingerprint=" << fingerprint << std::endl;
  }

}

#endif
