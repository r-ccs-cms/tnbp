/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/truncation.h
@brief truncation based on the messenger tensors
*/
#ifndef TNBP_PTNS_FUNCTION_TRUNCATION_H
#define TNBP_PTNS_FUNCTION_TRUNCATION_H

#include <type_traits>
#include "tnbp/framework/root.h"
#include "tnbp/framework/mpiutility.h"

namespace tnbp {

  /**
     Truncation based on the messenger tensors
   */
  template <typename TenT>
  void Truncation(context_handle_t<TenT> & ctx,
		  const std::vector<std::pair<int,int>> & I,
		  std::vector<TenT> & V,
		  const std::vector<int> & SiteIdx,
		  const std::map<int,int> & Site_To_MpiRank,
		  std::vector<TenT> & E,
		  const std::vector<int> & EdgeIdx,
		  MPI_Comm comm,
		  bond_dim_t<TenT> max_dim,
		  real_t<TenT> eps,
		  real_t<TenT> err,
		  std::vector<bond_dim_t<TenT>> & res_bond_dim,
		  std::vector<real_t<TenT>> & res_trunc_err,
		  std::vector<real_t<TenT>> * site_norms = nullptr) {
    
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using OrderT = typename tci::tensor_traits<TenT>::order_t;
    using SizeT = typename tci::tensor_traits<TenT>::ten_size_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tci::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tci::create_context(ctx_r);

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    size_t num_e = EdgeIdx.size();

    TenT Ra;
    TenT Sa;
    TenT Rb;
    TenT Sb;

    res_trunc_err.resize(num_e);
    res_bond_dim.resize(num_e);

    // Apply a two-leg update factor onto V[site] along its `bond_address` leg,
    // renormalize, and record the norm. `factor_first` selects the contraction
    // orientation: the site_a update keeps the factor on the right
    // (factor_first=false); the site_b update and the cross-rank receiver keep
    // it on the left (factor_first=true).
    auto apply_factor_to_site = [&](int site, int bond_address, TenT & factor,
				    bool factor_first) {
      auto it = std::find(SiteIdx.begin(),SiteIdx.end(),site);
      int site_address = std::distance(SiteIdx.begin(),it);
      OrderT ord = tci::order(ctx,V[site_address]);
      List<BondLabelT> IdxV(ord);
      List<BondLabelT> IdxC(ord);
      std::iota(IdxV.begin(),IdxV.end(),0);
      std::iota(IdxC.begin(),IdxC.end(),0);
      IdxV[bond_address] = static_cast<BondLabelT>(-1);
      List<BondLabelT> IdxF(2);
      if( factor_first ) {
	IdxF[0] = static_cast<BondLabelT>(bond_address);
	IdxF[1] = static_cast<BondLabelT>(-1);
	tci::contract(ctx,factor,IdxF,V[site_address],IdxV,V[site_address],IdxC);
      } else {
	IdxF[0] = static_cast<BondLabelT>(-1);
	IdxF[1] = static_cast<BondLabelT>(bond_address);
	tci::contract(ctx,V[site_address],IdxV,factor,IdxF,V[site_address],IdxC);
      }
      auto nrm = tci::normalize(ctx,V[site_address]);
      if (site_norms) site_norms->push_back(nrm);
    };

    for(int global_edge_address=0; global_edge_address < I.size(); global_edge_address++) {
      int site_a = I[global_edge_address].first;
      int site_b = I[global_edge_address].second;
      int mpi_rank_a = Site_To_MpiRank.at(site_a);
      int mpi_rank_b = Site_To_MpiRank.at(site_b);
      int mpi_type = 0;
      if( mpi_rank_a == mpi_rank ) {
	mpi_type += 2;
      }
      if( mpi_rank_b == mpi_rank ) {
	mpi_type += 1;
      }
      if( mpi_type > 0 ) {
	std::vector<int> bond_a = GetSurroundingBondIndex(site_a,I);
	std::vector<int> bond_b = GetSurroundingBondIndex(site_b,I);
	int target_edge = 0;
	int bond_address_a = 0;
	int bond_address_b = 0;
	int direction = 0;
	for(int k=0; k < bond_a.size(); k++) {
	  if( ( I[bond_a[k]].first == site_a ) &&
	      ( I[bond_a[k]].second == site_b ) ) {
	    target_edge = bond_a[k];
	    bond_address_a = k;
	    direction = 0;
	    break;
	  }
	  if( ( I[bond_a[k]].first == site_b ) &&
	      ( I[bond_a[k]].second == site_a ) ) {
	    target_edge = bond_a[k];
	    bond_address_a = k;
	    direction = 1;
	    break;
	  }
	}
	auto it_bond_address_b = std::find(bond_b.begin(),bond_b.end(),
					   target_edge);
	bond_address_b = std::distance(bond_b.begin(),it_bond_address_b);

	auto it_edge_address = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					 target_edge);
	int edge_address = std::distance(EdgeIdx.begin(),it_edge_address);

	// Cross-rank seam edges are truncated authoritatively on the site_a
	// owner (mpi_type==2): it runs the eigh/SVD pipeline once and ships the
	// site_b owner the factor to apply plus the spectrum and scalars. If both
	// incident ranks re-ran the pipeline on the (bit-identical) messenger they
	// could still disagree on the retained bond dimension, because eigh in
	// SquareRootAndInverse is not bit-reproducible across ranks on a
	// near-singular messenger and the amplified inverse-sqrt flips the
	// singular-value count near the cutoff. Receiving the result removes the
	// receiver's independent decomposition entirely.
	const bool cross_rank = (mpi_rank_a != mpi_rank_b);
	if( cross_rank && mpi_type == 1 ) {
	  TenT Mb;
	  MpiRecv(ctx,Mb,mpi_rank_a,comm);
	  RealTenT spec;
	  MpiRecv(ctx_r,spec,mpi_rank_a,comm);
	  RealT scalars[2] = {RealT(0),RealT(0)};
	  MPI_Recv(scalars,2,GetMpiType<RealT>(),mpi_rank_a,10,comm,
		   MPI_STATUS_IGNORE);
	  res_trunc_err[edge_address] = scalars[0];
	  RealT recv_norm_t = scalars[1];

	  // res_bond_dim is the retained singular-value count == spectrum length.
	  auto shape_spec = tci::shape(ctx_r,spec);
	  res_bond_dim[edge_address] = shape_spec[0];

	  // Rebuild E = diag(spectrum) with the same real/complex handling the
	  // authoritative rank uses for Z below.
	  TenT Z;
	  if constexpr (std::is_same_v<TenT,RealTenT>) {
	    Z = tci::copy(ctx_r,spec);
	  } else {
	    Z = tci::to_cplx(ctx_r,spec);
	  }
	  tci::diag(ctx,Z);
	  E[edge_address] = tci::copy(ctx,Z);
	  E[edge_address+num_e] = tci::copy(ctx,Z);

	  // Record norm_t (the authoritative rank's value) then norm_b, keeping
	  // the two-entries-per-seam order the consumer's cross-rank norm dedup
	  // relies on; that dedup discards this mpi_type==1 norm_t downstream.
	  if (site_norms) site_norms->push_back(recv_norm_t);
	  apply_factor_to_site(site_b,bond_address_b,Mb,/*factor_first=*/true);
	  continue; // seam handled from the authoritative rank's result
	}

	if( direction == 0 ) {
	  SquareRootAndInverse(ctx,E[edge_address],Ra,Sa,eps);
	  SquareRootAndInverse(ctx,E[edge_address+num_e],Rb,Sb,eps);
	} else {
	  SquareRootAndInverse(ctx,E[edge_address+num_e],Ra,Sa,eps);
	  SquareRootAndInverse(ctx,E[edge_address],Rb,Sb,eps);
	}

	TenT T;
	List<BondLabelT> IdxRa(2);
	List<BondLabelT> IdxRb(2);
	List<BondLabelT> IdxRR(2);
	IdxRa[0] = static_cast<BondLabelT>(-1);
	IdxRa[1] = static_cast<BondLabelT>(0);
	IdxRb[0] = static_cast<BondLabelT>(-1);
	IdxRb[1] = static_cast<BondLabelT>(1);
	IdxRR[0] = static_cast<BondLabelT>(0);
	IdxRR[1] = static_cast<BondLabelT>(1);
	tci::contract(ctx,Ra,IdxRa,Rb,IdxRb,T,IdxRR);
	auto norm_t = tci::normalize(ctx,T);
	if (site_norms) site_norms->push_back(norm_t);

	TenT X;
	TenT Y;
	RealTenT S;
	OrderT num_rows = 1;
	RealT trunc_err;
	BondDimT chi_min = 1;
	BondDimT chi_max = max_dim;
	tci::trunc_svd(ctx,T,num_rows,X,S,Y,
		       trunc_err,chi_min,chi_max,err,eps);
	
	auto shape_s = tci::shape(ctx_r,S);
	res_bond_dim[edge_address] = shape_s[0];
	res_trunc_err[edge_address] = trunc_err;
	// Snapshot the pristine spectrum before the in-place sqrt below overwrites
	// S, but only on the cross-rank sender: it ships the spectrum so the
	// receiver can rebuild diag(spectrum) for its E slots. Intra-rank and
	// single-rank edges never send, so they skip this copy in the hot loop.
	RealTenT spec;
	if( cross_rank && mpi_type == 2 ) spec = tci::copy(ctx_r,S);
	TenT Z;
	TenT P;
	if constexpr (std::is_same_v<TenT,RealTenT>) {
	  Z = tci::copy(ctx_r,S);
	  tci::for_each(ctx_r,S,[](auto & elem){ elem = std::sqrt(elem); });
	  tci::move(ctx_r,S,P);
	} else {
	  Z = tci::to_cplx(ctx_r,S);
	  tci::for_each(ctx_r,S,[](auto & elem){ elem = std::sqrt(elem); });
	  P = tci::to_cplx(ctx_r,S);
	}
	tci::diag(ctx,Z);
	tci::diag(ctx,P);
	E[edge_address] = tci::copy(ctx,Z);
	E[edge_address+num_e] = tci::copy(ctx,Z);

	List<BondLabelT> IdxU(2);
	List<BondLabelT> IdxS(2);
	List<BondLabelT> IdxP(2);
	List<BondLabelT> IdxT(2);

	if( mpi_type == 2 || mpi_type == 3 ) {
	  IdxS[0] = static_cast<BondLabelT>(-1);
	  IdxS[1] = static_cast<BondLabelT>(0);
	  IdxU[0] = static_cast<BondLabelT>(-1);
	  IdxU[1] = static_cast<BondLabelT>(1);
	  IdxT[0] = static_cast<BondLabelT>(0);
	  IdxT[1] = static_cast<BondLabelT>(1);
	  tci::contract(ctx,Sa,IdxS,X,IdxU,T,IdxT);
	  IdxT[0] = static_cast<BondLabelT>(0);
	  IdxT[1] = static_cast<BondLabelT>(-1);
	  IdxP[0] = static_cast<BondLabelT>(-1);
	  IdxP[1] = static_cast<BondLabelT>(1);
	  IdxU[0] = static_cast<BondLabelT>(0);
	  IdxU[1] = static_cast<BondLabelT>(1);
	  tci::contract(ctx,T,IdxT,P,IdxP,T,IdxU);
	  apply_factor_to_site(site_a,bond_address_a,T,/*factor_first=*/false);
	}

	if( mpi_type == 1 || mpi_type == 3 ) {
	  IdxU[0] = static_cast<BondLabelT>(0);
	  IdxU[1] = static_cast<BondLabelT>(-1);
	  IdxS[0] = static_cast<BondLabelT>(-1);
	  IdxS[1] = static_cast<BondLabelT>(1);
	  IdxT[0] = static_cast<BondLabelT>(0);
	  IdxT[1] = static_cast<BondLabelT>(1);
	  tci::contract(ctx,Y,IdxU,Sb,IdxS,T,IdxT);
	  IdxP[0] = static_cast<BondLabelT>(0);
	  IdxP[1] = static_cast<BondLabelT>(-1);
	  IdxT[0] = static_cast<BondLabelT>(-1);
	  IdxT[1] = static_cast<BondLabelT>(1);
	  IdxU[0] = static_cast<BondLabelT>(0);
	  IdxU[1] = static_cast<BondLabelT>(1);
	  tci::contract(ctx,P,IdxP,T,IdxT,T,IdxU);
	  apply_factor_to_site(site_b,bond_address_b,T,/*factor_first=*/true);
	}

	if( cross_rank && mpi_type == 2 ) {
	  // Build the site_b owner's update factor M_b = P·(Y·Sb) with the same
	  // contractions the mpi_type==1 block runs, then ship it with the
	  // pristine spectrum and the (trunc_err, norm_t) scalars so the receiver
	  // applies an identical, authoritative result. IdxU/IdxS/IdxP/IdxT are
	  // the reusable label buffers declared before the mpi_type blocks.
	  TenT Tb;
	  IdxU[0] = static_cast<BondLabelT>(0);
	  IdxU[1] = static_cast<BondLabelT>(-1);
	  IdxS[0] = static_cast<BondLabelT>(-1);
	  IdxS[1] = static_cast<BondLabelT>(1);
	  IdxT[0] = static_cast<BondLabelT>(0);
	  IdxT[1] = static_cast<BondLabelT>(1);
	  tci::contract(ctx,Y,IdxU,Sb,IdxS,Tb,IdxT);
	  TenT Mb;
	  IdxP[0] = static_cast<BondLabelT>(0);
	  IdxP[1] = static_cast<BondLabelT>(-1);
	  IdxT[0] = static_cast<BondLabelT>(-1);
	  IdxT[1] = static_cast<BondLabelT>(1);
	  IdxU[0] = static_cast<BondLabelT>(0);
	  IdxU[1] = static_cast<BondLabelT>(1);
	  tci::contract(ctx,P,IdxP,Tb,IdxT,Mb,IdxU);
	  MpiSend(ctx,Mb,mpi_rank_b,comm);
	  MpiSend(ctx_r,spec,mpi_rank_b,comm);
	  RealT scalars[2] = { res_trunc_err[edge_address], norm_t };
	  MPI_Send(scalars,2,GetMpiType<RealT>(),mpi_rank_b,10,comm);
	}
      }
    }

    
  }

  
  
}

#endif
