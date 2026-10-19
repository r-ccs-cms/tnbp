/// This file is a part of r-ccs-cms/tnbp
/**
   @file tnbp/ptns/function/bp.h
   @brief Belief propagation methods
*/
#ifndef TNBP_PTNS_FUNCTION_BP_H
#define TNBP_PTNS_FUNCTION_BP_H

#include "tnbp/ptns/tps.h"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/graph.h"
#include "tnbp/framework/mpiutility.h"
#include "tnbp/ptns/function/message_gauge.h"

namespace tnbp {

/** Build a local boundary-edge schedule from disjoint rank-pair rounds.
 * All ranks derive identical rounds from the shared graph/ownership. Each
 * round greedily selects rank pairs without reusing a rank, and exchanges all
 * edges of each selected pair in global edge order. Only this rank's edges
 * are retained. Rebuild when graph or ownership changes, not bond dimensions.
 */
inline std::vector<int> MakeBPCommunicationOrder(
    const std::vector<std::pair<int,int>>& I,
    const std::map<int,int>& owners, int rank) {
  std::map<std::pair<int,int>,std::vector<int>> pending;
  for(std::size_t edge=0; edge<I.size(); ++edge) {
    const int a=owners.at(I[edge].first), b=owners.at(I[edge].second);
    if(a!=b) pending[{std::min(a,b),std::max(a,b)}].push_back(static_cast<int>(edge));
  }
  std::vector<int> result;
  while(!pending.empty()) {
    std::set<int> used;
    for(auto it=pending.begin(); it!=pending.end();) {
      const auto pair=it->first;
      if(used.count(pair.first) || used.count(pair.second)) {
        ++it;
        continue;
      }
      used.insert(pair.first);
      used.insert(pair.second);
      if(rank==pair.first || rank==pair.second)
        result.insert(result.end(),it->second.begin(),it->second.end());
      it=pending.erase(it);
    }
  }
  return result;
}

/** Compute all locally sourced messages from immutable E into F.
 * F uses the same two-direction layout as E; remote-source slots are filled
 * by ExchangeBPMessages. Intermediate normalization follows the legacy update.
 */
template<typename TenT>
void UpdateBPMessages(context_handle_t<TenT>& ctx,
    const std::vector<std::pair<int,int>>& I, const std::vector<TenT>& V,
    const std::vector<int>& SiteIdx, const std::vector<TenT>& E,
    const std::vector<int>& EdgeIdx, MPI_Comm comm, std::vector<TenT>& F) {
  using BondLabelT=typename tcapi::tensor_traits<TenT>::bond_label_t;
  int rank; MPI_Comm_rank(comm,&rank);
  static bool logged=false;
  if(!logged) {
    if(rank==0) std::clog << "BP_ALGORITHM " << BPAlgorithm << " gauge=" << BPGauge
                         << " residual=" << BPResidual << std::endl;
    logged=true;
  }
  const std::size_t size_e=EdgeIdx.size();
  F.resize(2*size_e);
  for(std::size_t address=0; address<SiteIdx.size(); ++address) {
    const int site=SiteIdx[address];
    const auto bonds=GetSurroundingBondIndex(site,I);
    for(std::size_t target=0; target<bonds.size(); ++target) {
      const int edge=bonds[target];
      TenT message;
      TenT T = tcapi::copy(ctx,V[address]);
      auto order_v = tcapi::order(ctx,T);
      List<BondLabelT> idx_v(order_v);
      List<BondLabelT> idx_e(2);
      List<BondLabelT> idx_result(order_v);

      for(int k=0; k < bonds.size(); k++) {
        auto it_message = std::find(EdgeIdx.begin(),EdgeIdx.end(),
                               bonds[k]);
        int message_slot = std::distance(EdgeIdx.begin(),it_message);
        if( k != target ) {
          std::iota(idx_v.begin(),idx_v.end(),0);
          std::iota(idx_result.begin(),idx_result.end(),0);
          idx_v[k] = static_cast<BondLabelT>(-1);
          idx_e[0] = static_cast<BondLabelT>(-1);
          idx_e[1] = static_cast<BondLabelT>(k);
          if( I[bonds[k]].first == site ) {
            tcapi::contract(ctx,T,idx_v,E[message_slot+size_e],idx_e,T,idx_result);
          } else {
            tcapi::contract(ctx,T,idx_v,E[message_slot],idx_e,T,idx_result);
          }
        }
        tcapi::normalize(ctx,T);
      }
      message = tcapi::copy(ctx,V[address]);
      tcapi::cplx_conj(ctx,message);

      List<BondLabelT> idx_t(order_v);
      List<BondLabelT> idx_message(2);
      std::iota(idx_v.begin(),idx_v.end(),-order_v);
      std::iota(idx_t.begin(),idx_t.end(),-order_v);
      std::iota(idx_message.begin(),idx_message.end(),0);
      idx_t[target] = static_cast<BondLabelT>(0);
      idx_v[target] = static_cast<BondLabelT>(1);
      tcapi::contract(ctx,T,idx_t,message,idx_v,message,idx_message);
      NormalizeMessage(ctx,message,comm,"BP edge="+std::to_string(edge)+" source="+std::to_string(site));
      const auto pos=std::find(EdgeIdx.begin(),EdgeIdx.end(),edge);
      const std::size_t slot=std::distance(EdgeIdx.begin(),pos)
          +(I[edge].second==site ? size_e : 0);
      F[slot]=tcapi::move(ctx,message);
    }
  }
}

/** Exchange only boundary messages after local computation has completed.
 * order must come from MakeBPCommunicationOrder: disjoint rank pairs within
 * each shared round, with matching edge order at both endpoints. The first
 * endpoint sends then receives; the second receives then sends. Blocking
 * dependencies cannot cycle; independent pairs may communicate concurrently.
 * No round barrier is needed (rounds may overlap in time).
 */
template<typename TenT>
void ExchangeBPMessages(context_handle_t<TenT>& ctx,
    const std::vector<std::pair<int,int>>& I,
    const std::map<int,int>& owners, const std::vector<int>& EdgeIdx,
    const std::vector<int>& order, MPI_Comm comm, std::vector<TenT>& F) {
  int rank; MPI_Comm_rank(comm,&rank);
  const std::size_t n=EdgeIdx.size();
  for(const int edge:order) {
    const int a=owners.at(I[edge].first), b=owners.at(I[edge].second);
    const std::size_t slot=std::distance(EdgeIdx.begin(),
        std::find(EdgeIdx.begin(),EdgeIdx.end(),edge));
    if(rank==a) {
      MpiSend(ctx,F[slot],b,comm);
      MpiRecv(ctx,F[slot+n],b,comm);
    } else {
      MpiRecv(ctx,F[slot],a,comm);
      MpiSend(ctx,F[slot+n],a,comm);
    }
  }
}

/** Compatibility entry point for subset-edge updates.
 * Computes local predictions into a temporary buffer, but exchanges and commits
 * only J. Repeated full sweeps should use UpdateBPMessages/ExchangeBPMessages
 * directly to avoid recomputing predictions for each legacy layer.
 */
  template <typename TenT>
  void BeliefPropagation(context_handle_t<TenT> & ctx,
			 const std::vector<std::pair<int,int>> & I,
			 const std::vector<TenT> & V,
			 const std::vector<int> & SiteIdx,
			 const std::map<int,int> & Site_To_MpiRank,
			 const std::vector<std::pair<int,int>> & J,
			 const std::vector<TenT> & E,
			 const std::vector<int> & EdgeIdx,
			 MPI_Comm comm,
			 std::vector<TenT> & F) {

    std::vector<TenT> all;
    UpdateBPMessages(ctx,I,V,SiteIdx,E,EdgeIdx,comm,all);
    int rank; MPI_Comm_rank(comm,&rank);
    auto order=MakeBPCommunicationOrder(I,Site_To_MpiRank,rank);
    auto selected=[&](int edge) {
      return std::any_of(J.begin(),J.end(),[&](const auto& pair) {
        return pair==I[edge] || (pair.first==I[edge].second && pair.second==I[edge].first);
      });
    };
    order.erase(std::remove_if(order.begin(),order.end(),
        [&](int edge){return !selected(edge);}),order.end());
    ExchangeBPMessages(ctx,I,Site_To_MpiRank,EdgeIdx,order,comm,all);
    F.resize(2*EdgeIdx.size());
    for(std::size_t j=0;j<EdgeIdx.size();++j) if(selected(EdgeIdx[j])) {
      F[j]=tcapi::move(ctx,all[j]);
      F[j+EdgeIdx.size()]=tcapi::move(ctx,all[j+EdgeIdx.size()]);
    }
  }

  /**
     Legacy fixed-point check: reconstruct predictions from current messages.
     Kept for compatibility and diagnostic comparisons. New BP loops use
     BeliefPropagationDifference on the actual previous/updated pair.
   */
  template <typename TenT>
  void BeliefPropagationCondition(context_handle_t<TenT> & ctx,
				  const std::vector<std::pair<int,int>> & I,
				  const std::vector<TenT> & V,
				  const std::vector<int> & SiteIdx,
				  const std::map<int,int> & Site_To_MpiRank,
				  const std::vector<TenT> & E,
				  const std::vector<int> & EdgeIdx,
				  MPI_Comm comm,
				  real_t<TenT> & result) {

    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using RealT = typename tcapi::tensor_traits<TenT>::real_t;
    using RealTenT = typename tcapi::tensor_traits<TenT>::real_ten_t;
    using OrderT = typename tcapi::tensor_traits<TenT>::order_t;
    using BondLabelT  = typename tcapi::tensor_traits<TenT>::bond_label_t;

    std::size_t num_v = SiteIdx.size();
    std::size_t num_e = EdgeIdx.size();
    std::size_t num_total_edges = I.size();
    RealT result_volume = 1.0/(2.0*num_total_edges);
    RealT result_rank = 0.0;
    result = 0.0;

    for(std::size_t address=0; address < SiteIdx.size(); address++) {
      std::vector<int> BondIdx = GetSurroundingBondIndex(SiteIdx[address],I);
      for(std::size_t k=0; k < BondIdx.size(); k++) {
	TenT T = tcapi::copy(ctx,V[address]);
	auto OrderV = tcapi::order(ctx,T);
	for(std::size_t l=0; l < BondIdx.size(); l++) {
	  if( l != k ) {
	    auto itEdgeAddress = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					   BondIdx[l]);
	    int EdgeAddress = std::distance(EdgeIdx.begin(),itEdgeAddress);
	    if( I[BondIdx[l]].first == SiteIdx[address] ) {
	      EdgeAddress += num_e;
	    }
	    const TenT & F = E[EdgeAddress];
	    List<BondLabelT> IdxE(2);
	    List<BondLabelT> IdxV(OrderV);
	    List<BondLabelT> IdxV_res(OrderV);
	    std::iota(IdxV.begin(),IdxV.end(),0);
	    std::iota(IdxV_res.begin(),IdxV_res.end(),0);
	    IdxV[l] = static_cast<BondLabelT>(-1);
	    IdxE[0] = static_cast<BondLabelT>(-1);
	    IdxE[1] = static_cast<BondLabelT>(l);
	    tcapi::contract(ctx,T,IdxV,F,IdxE,T,IdxV_res);
	    auto norm_v = tcapi::normalize(ctx,T);
	  }
	} // end for(std::size_t l=0; l < BondIdx.size(); l++)

#ifdef BP_DEBUG
	std::cout << " belief propagation: Before physical contruction " << std::endl;
#endif

	TenT C = tcapi::copy(ctx,V[address]);
	tcapi::cplx_conj(ctx,C);
	List<BondLabelT> IdxV(OrderV);
	List<BondLabelT> IdxC(OrderV);
	List<BondLabelT> IdxN(2);
	std::iota(IdxV.begin(),IdxV.end(),-OrderV);
	std::iota(IdxC.begin(),IdxC.end(),-OrderV);
	IdxV[k] = static_cast<BondLabelT>(0);
	IdxC[k] = static_cast<BondLabelT>(1);
	IdxN[0] = static_cast<BondLabelT>(0);
	IdxN[1] = static_cast<BondLabelT>(1);
	tcapi::contract(ctx,T,IdxV,C,IdxC,T,IdxN);
	auto norm_v = tcapi::normalize(ctx,T);

	auto itEdgeAddress = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				       BondIdx[k]);
	int EdgeAddress = std::distance(EdgeIdx.begin(),itEdgeAddress);
	if( I[BondIdx[k]].second == SiteIdx[address] ) {
	  EdgeAddress += num_e;
	}
	TenT F = tcapi::copy(ctx,E[EdgeAddress]);
	auto norm_f = tcapi::normalize(ctx,F);

#ifdef BP_DEBUG
	std::cout << " belief propagation: Before linear_combine " << std::endl;
#endif

	FixMessageTracePhase(ctx,T,"BP prediction edge="+std::to_string(BondIdx[k]),comm);
        FixMessageTracePhase(ctx,F,"BP reference edge="+std::to_string(BondIdx[k]),comm);
	TenT D = tcapi::linear_combine<TenT>(ctx,{std::cref(T),std::cref(F)},{ElemT(1.0),ElemT(-1.0)});
	RealT NormD = tcapi::norm(ctx,D);
	result_rank += NormD * result_volume;

      } // end for(std::size_t k=0; k < BondIdx.size(); k++)
    } // end for(int address=0; address < SiteIdx.size(); address++)
    MPI_Allreduce(&result_rank,&result,1,GetMpiType<RealT>(),MPI_SUM,comm);
  }


  /**
   * @brief Mean directed-edge distance between successive BP messages.
   *
   * Both inputs must satisfy the messenger invariant: unit Frobenius norm
   * and positive-real trace. Only the source site's owner contributes each
   * directed edge, so inter-rank replica slots are not counted twice. The
   * scalar sum is reduced across comm and divided by 2*I.size(). No site
   * tensors or new BP predictions are computed here.
   *
   * Call after the entire synchronous update has filled updated, before
   * replacing previous. On convergence callers still commit updated. This
   * stopping rule differs from the legacy BeliefPropagationCondition, which
   * reconstructs another prediction from the already updated messages.
   */
  template<class TenT>
  void BeliefPropagationDifference(context_handle_t<TenT>& ctx,
      const std::vector<std::pair<int,int>>& I,
      const std::map<int,int>& Site_To_MpiRank,
      const std::vector<TenT>& previous,
      const std::vector<TenT>& updated,
      const std::vector<int>& EdgeIdx,
      MPI_Comm comm, real_t<TenT>& result) {
    using ElemT = elem_t<TenT>;
    using RealT = real_t<TenT>;
    const std::size_t n = EdgeIdx.size();
    if(previous.size()!=2*n || updated.size()!=2*n)
      throw std::invalid_argument("BP difference: message count must equal twice the local edge count");
    int rank; MPI_Comm_rank(comm,&rank);
    static bool logged=false;
    if(!logged) {
      if(rank==0) std::clog << "BP_CONVERGENCE successive-normalized-message-difference-v1" << std::endl;
      logged=true;
    }
    RealT local=0;
    const RealT volume=I.empty() ? RealT(0) : RealT(1)/(RealT(2)*RealT(I.size()));
    for(std::size_t j=0;j<n;++j) {
      const auto& edge=I.at(EdgeIdx[j]);
      for(int direction=0;direction<2;++direction) {
        const int source=direction ? edge.second : edge.first;
        if(Site_To_MpiRank.at(source)!=rank) continue;
        const std::size_t slot=j+direction*n;
        TenT delta=tcapi::linear_combine<TenT>(ctx,
            {std::cref(updated[slot]),std::cref(previous[slot])},{ElemT(1),ElemT(-1)});
        local+=tcapi::norm(ctx,delta)*volume;
      }
    }
    MPI_Allreduce(&local,&result,1,GetMpiType<RealT>(),MPI_SUM,comm);
  }

  /** Build the boundary schedule for a TensorProductState. */
  template <typename TenT>
  std::vector<int> MakeBPCommunicationOrder(const TensorProductState<TenT>& W) {
    return MakeBPCommunicationOrder(W.I_,W.Site_To_MpiRank_,W.mpi_rank_);
  }

  /** Compute local messages using the state-owned communicator. */
  template <typename TenT>
  void UpdateBPMessages(context_handle_t<TenT>& ctx,
                        const TensorProductState<TenT>& W,
                        const std::vector<TenT>& E,
                        const std::vector<int>& EdgeIdx,
                        std::vector<TenT>& F) {
    UpdateBPMessages(ctx,W.I_,W.V_,W.SiteIdx_,E,EdgeIdx,W.comm_,F);
  }

  /** Complete boundary messages using the state-owned ownership map. */
  template <typename TenT>
  void ExchangeBPMessages(context_handle_t<TenT>& ctx,
                          const TensorProductState<TenT>& W,
                          const std::vector<int>& EdgeIdx,
                          const std::vector<int>& order,
                          std::vector<TenT>& F) {
    ExchangeBPMessages(ctx,W.I_,W.Site_To_MpiRank_,EdgeIdx,order,W.comm_,F);
  }

  /** Legacy subset-update wrapper retained for external applications.
   * New callers should compute local messages and exchange boundaries explicitly.
   */
  template <typename TenT>
  void BeliefPropagation(context_handle_t<TenT> & ctx,
			 const TensorProductState<TenT> & W,
			 const std::vector<std::pair<int,int>> & J,
			 const std::vector<TenT> & E,
			 const std::vector<int> & EdgeIdx,
			 std::vector<TenT> & F) {
    BeliefPropagation(ctx,W.I_,W.V_,W.SiteIdx_,W.Site_To_MpiRank_,
			   J,E,EdgeIdx,W.comm_,F);
  }

  /**
     Legacy fixed-point check wrapper for TensorProductState.
     Reconstructs predictions; it is not the successive-message difference.
   */
  template <typename TenT>
  void BeliefPropagationCondition(context_handle_t<TenT> & ctx,
				  const TensorProductState<TenT> & W,
				  const std::vector<TenT> & E,
				  const std::vector<int> & EdgeIdx,
				  real_t<TenT> & result) {
    BeliefPropagationCondition(ctx,W.I_,W.V_,W.SiteIdx_,W.Site_To_MpiRank_,
			       E,EdgeIdx,W.comm_,result);
  }


}

#endif

