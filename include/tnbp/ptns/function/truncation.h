/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/truncation.h
@brief truncation based on the messenger tensors
*/
#ifndef TNBP_PTNS_FUNCTION_TRUNCATION_H
#define TNBP_PTNS_FUNCTION_TRUNCATION_H

#include "tnbp/framework/root.h"

namespace tnbp {

  /**
     Truncation based on the messenger tensors
   */
  template <typename TenT>
  void Truncation(context_handle_t<TenT> & ctx,
		  const std::vector<std::pair<int,int>> & I,
		  std::vector<TenT> & V,
		  const std::vector<int> & SiteIdx,
		  const std::vector<int> & Site_To_MpiRank,
		  std::vector<TenT> & E,
		  const std::vector<int> & EdgeIdx,
		  MPI_Comm comm,
		  bond_dim_t<TenT> max_dim,
		  real_t<TenT> eps,
		  real_t<TenT> err) {
    
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
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

    for(int edge_address=0; edge_address < I.size(); edge_address++) {
      int site_a = I[edge_address].first;
      int site_b = I[edge_address].second;
      int mpi_rank_a = Site_To_MpiRank[site_a];
      int mpi_rank_b = Site_To_MpiRank[site_b];
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

	int site_address_a;
	int site_address_b;
	auto it_edge_address = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					 target_edge);
	int edge_address = std::distance(EdgeIdx.begin(),it_edge_address);

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
	tci::contract(ctx,Ra,IdxRa,Rb,IdxRb,T,IdxRR);
	auto norm_t = tci::normalize(ctx,T);

	TenT X;
	TenT Y;
	RealTenT S;
	RankT num_rows = 1;
	RealT trunc_err;
	BondDimT chi_min = 1;
	BondDimT chi_max = max_dim;
	tci::trunc_svd(ctx,T,num_rows,X,S,Y,
		       trunc_err,chi_min,chi_max,err,eps);

	TenT Z;
	tci::convert(ctx_r,S,ctx,Z);
	tci::for_each(ctx_r,S,[](RealT & elem){ elem = std::sqrt(elem); });
	TenT P;
	tci::convert(ctx_r,S,ctx,P);
	tci::diag(ctx,Z);
	tci::diag(ctx,P);
	tci::copy(ctx,Z,E[edge_address]);

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
	  auto it_site_address_a = std::find(SiteIdx.begin(),SiteIdx.end(),
					site_a);
	  site_address_a = std::distance(SiteIdx.begin(),it_site_address_a);

	  RankT rank_a = tci::rank(ctx,V[site_address_a]);
	  List<BondLabelT> IdxA(rank_a);
	  List<BondLabelT> IdxC(rank_a);
	  std::iota(IdxA.begin(),IdxA.end(),0);
	  std::iota(IdxC.begin(),IdxC.end(),0);
	  IdxA[bond_address_a] = static_cast<BondLabelT>(-1);
	  IdxT[0] = static_cast<BondLabelT>(-1);
	  IdxT[1] = static_cast<BondLabelT>(bond_address_a);
	  tci::contract(ctx,V[site_address_a],IdxA,T,IdxT,
			V[site_address_a],IdxC);
	  auto norm_a = tci::normalize(ctx,V[site_address_a]);
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
	  auto it_site_address_b = std::find(SiteIdx.begin(),SiteIdx.end(),
					     site_b);
	  site_address_b = std::distance(SiteIdx.begin(),it_site_address_b);
	  RankT rank_b = tci::rank(ctx,V[site_address_b]);
	  List<BondLabelT> IdxB(rank_b);
	  List<BondLabelT> IdxC(rank_b);
	  std::iota(IdxB.begin(),IdxB.end(),0);
	  std::iota(IdxC.begin(),IdxC.end(),0);
	  IdxB[bond_address_b] = static_cast<BondLabelT>(-1);
	  IdxT[0] = static_cast<BondLabelT>(bond_address_b);
	  IdxT[1] = static_cast<BondLabelT>(-1);
	  tci::contract(ctx,T,IdxT,V[site_address_b],IdxB,
			V[site_address_b],IdxC);
	  auto norm_b = tci::normalize(ctx,V[site_address_b]);
	}
      }
    }

    
  }
  
}

#endif
