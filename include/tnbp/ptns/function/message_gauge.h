#pragma once

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/mpiutility.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace tnbp {
inline constexpr const char* BPAlgorithm = "trace-phase-gauge-v1";
inline constexpr const char* BPGauge = "positive-real-trace;threshold=64*epsilon*frobenius-norm";
inline constexpr const char* BPResidual = "mean-directed-edge-frobenius-distance-of-unit-norm-trace-gauged-messages";

// Select a scalar-phase representative only; this is not a PSD projection.
// A null communicator allows callers to catch domain_error. Distributed BP
// aborts its communicator on invalid messages, avoiding stranded peer receives.
template<class TenT>
void FixMessageTracePhase(context_handle_t<TenT>& ctx, TenT& message,
                          const std::string& location = "message",
                          MPI_Comm comm = MPI_COMM_NULL) {
  try {
  using Elem = elem_t<TenT>;
  using Real = real_t<TenT>;
  const auto shape = tcapi::shape(ctx,message);
  if (shape.size()!=2 || shape[0]!=shape[1])
    throw std::invalid_argument("trace-phase gauge requires a square rank-2 message");
  const Real norm = tcapi::norm(ctx,message);
  TenT traced;
  tcapi::trace(ctx,message,{{0,1}},traced);
  // Empty coordinates address the logical rank-0 value, including lazy scale.
  const Elem value = tcapi::get_elem(ctx,traced,{});
  const Real magnitude = std::abs(value);
  const Real threshold = Real(64)*std::numeric_limits<Real>::epsilon()*norm;
  if (!std::isfinite(norm) || !(norm>Real(0)) ||
      !std::isfinite(magnitude) || magnitude<=threshold) {
    std::ostringstream diagnostic;
    diagnostic << std::setprecision(std::numeric_limits<Real>::max_digits10)
               << BPAlgorithm << ": undefined trace phase at " << location
               << " norm=" << norm << " trace=" << value << " threshold=" << threshold;
    throw std::domain_error(diagnostic.str());
  }
  Elem factor;
  if constexpr (std::is_arithmetic_v<Elem>) factor=value>0 ? Elem(1) : Elem(-1);
  else factor=std::conj(value)/magnitude;
  tcapi::scale(ctx,message,factor);
  } catch (const std::exception& error) {
    if (comm!=MPI_COMM_NULL) {
      int rank=0,size=1;
      MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
      std::cerr << BPAlgorithm << " rank=" << rank << " at " << location
                << ": " << error.what() << std::endl;
      if(size>1) MPI_Abort(comm,1);
    }
    throw;
  }
}

template<class TenT>
void FixMessageTraceGauges(context_handle_t<TenT>& ctx, std::vector<TenT>& messages,
                           const std::vector<int>& edges) {
  if (messages.size()!=2*edges.size())
    throw std::invalid_argument("trace-phase gauge: message/edge count mismatch");
  for (std::size_t m=0;m<messages.size();++m)
    FixMessageTracePhase(ctx,messages[m],"checkpoint edge="+std::to_string(edges[m%edges.size()])+
                         " direction="+std::to_string(m/edges.size()));
}
}
