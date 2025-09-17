/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function.measurement.h
@brief Functions for measurement 
 */
#ifndef TNBP_PTNS_FUNCTION_MEASURE_H
#define TNBP_PTNS_FUNCTION_MEASURE_H

namespace tnbp {

  /**
     Function to get operator expectation value defined on the single site
   */
  template <typename TenT>
  std::vector<elem_t<TenT>> Measure(context_handle_t<TenT> & ctx,
			    const std::vector<std::pair<int,int>> & I,
			    const std::vector<TenT> & V,
			    const std::vector<int> & SiteIdx,
			    const std::vector<int> & Site_To_MpiRank,
			    const std::vector<TenT> & E,
			    const std::vector<int> EdgeIdx,
			    MPI_Comm comm,
			    const std::vector<int> & Site,
			    const std::vector<TenT> & O) {

    using ElemT = typename tci:;tensor_traits<TenT>::elem_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t size_e = EdgeIdx.size();

    size_t num_meas = 0;
    for(int site_address=0; site_address < Site.size(); site_address++) {
      if( mpi_rank == Site_To_MpiRank[Site[site_address]] ) {
	num_meas++;
      }
    }

    std::vector<elem_t<TenT>> result(num_meas,elem_t<TenT>(0.0));

    for(int site_address=0; site_address < Site; site_address++) {
      if( mpi_rank == Site_To_MpiRank[Site[site_address]] ) {
	auto it_siteidx_address = std::find(SiteIdx.begin(),SiteIdx.end(),
					    Site[site_address]);
	auto siteidx_address = std::distance(SiteIdx.begin(),it_siteidx_address);
	TenT W;
	tci::copy(ctx,V[siteidx_address],W);
	TenT Wdag;
	tci::copy(ctx,V[siteidx_address],Wdag);
	tci::cplx_conj(ctx,Wdag);
	RankT rank_w = tci::rank(ctx,W);
	std::vector<BondLabelT> IdxW(rank_w);
	std::vector<BondLabelT> IdxC(rank_w);
	std::vector<int> BondIdx = GetSurroundingBondIndex(Site[site_address],I);
	for(size_t m=0; m < BondIdx.size(); m++) {
	  auto it_edgeidx_address = std::find(EdgeIdx.begin(),EdgeIdx.end(),BondIdx[m]);
	  auto edgeidx_address = std::distance(EdgeIdx.begin(),it_edgeidx_address);
	  TenT F;
	  tci::copy(ctx,E[edgeidx_address+size_e],F);
	  std::vector<BondLabelT> IdxE(2);
	  std::vector<BondLabelT> IdxF(2);
	  std::vector<BondLabelT> IdxR(2);
	  IdxE[0] = 0;
	  IdxE[1] = -1;
	  IdxF[0] = -1;
	  IdxF[1] = 1;
	  IdxR[0] = 0;
	  IdxR[1] = 1;
	  tci::contract(ctx,F,IdxF,E[edgeidx_address+size_e],IdxE,F,IdxR);
	  std::iota(IdxW.begin(),IdxW.end(),0);
	  std::iota(IdxC.begin(),IdxC.end(),0);
	  IdxW[m] = -1;
	  IdxF[0] = -1;
	  IdxF[1] = m;
	  tci::contract(ctx,W,IdxW,F,IdxF,W,IdxC);
	}
	TenT R;
	std::vector<BondLabelT> IdxO(2);
	std::iota(IdxW.begin(),IdxW.end(),0);
	std::iota(IdxC.begin(),IdxC.end(),0);
	IdxW[rank_w-1] = -1;
	IdxO[0] = rank_w-1;
	IdxO[1] = -1;
	tci::contract(ctx,W,IdxW,O[i],IdxO,&R,IdxC);
	std::iota(IdxW.begin(),IdxW.end(),-rank_w);
	std::iota(IdxC.begin(),IdxC.end(),-rank_w);
	std::vector<BondLabelT> IdxM;
	tci::contract(ctx,W,IdxW,C,IdxC,W,IdxM);
	tci::contract(ctx,R,IdxW,C,IdxC,R,IdxM);
	RankT rank_r = tci::rank(ctx,R);
	CoorsT coors(rank_r,0);
	ElemT res_r = tci::get_elem(ctx,R,coors);
	ElemT res_w = tci::get_elem(ctx,W,coors);
	results[site_address] = res_r / res_w;
      }
    }
  }

  /**
     Function for the measurements of nearest-neighbor two-site local operators
   */
  template <typename TenT>
  std::vector<elem_t<TenT>> Measure(context_handle_t<TenT> & ctx,
			    const std::vector<std::pair<int,int>> & I,
			    const std::vector<TenT> & V,
			    const std::vector<int> & SiteIdx,
			    const std::vector<int> & Site_To_MpiRank,
			    const std::vector<TenT> & E,
			    const std::vector<int> EdgeIdx,
			    MPI_Comm comm,
			    const std::vector<std::pair<int,int>> & Edge,
			    const std::vector<TenT> & O) {

    using ElemT = typename tci:;tensor_traits<TenT>::elem_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t size_e = EdgeIdx.size();

    for(int m=0; m < Edge.size(); m++) {
      int site_a = Edge[m].first;
      int site_b = Edge[m].second;
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
	std::vector<int> bond_idx_a = GetSurroundingBondIndex(site_a,I);
	std::vector<int> bond_idx_b = GetSurroundingBondIndex(site_b,I);
	int target_edge_address = 0;
	int target_bond_address_a = 0;
	int target_bond_address_b = 0;
	for(int k=0; k < bond_idx_a.size(); k++) {
	  if( ( I[bond_idx_a[k]].first == site_a ) &&
	      ( I[bond_idx_a[k]].second == site_b ) ) {
	    target_edge_address = bond_idx_a[k];
	    target_bond_address_a = k;
	  }
	  if( ( I[bond_idx_a[k]].first == site_b ) &&
	      ( I[bond_idx_a[k]].second == site_a ) ) {
	    target_edge_address = bond_idx_a[k];
	    target_bond_address_a = k;
	  }
	}
	auto it_target_bond_address_b = std::find(bond_idx_b.begin(),bond_idx_b.end(),
						  target_edge_address);
	target_bond_address_b = std::distance(bond_idx_b.egin(),
					      it_target_bond_address_b);
	int site_address_a = 0;
	int site_address_b = 0;
	TenT A;
	TenT B;
	TenT AdagA;
	TenT BdagB;
	std::vector<BondLabelT> IdxE(2);

	if( mpi_type == 3 || mpi_type == 2 ) {
	  auot it_site_address_a = std::find(SiteIdx.begin(),SiteIdx.end(),site_a);
	  site_address_a = std::distance(SiteIdx.begin(),it_site_address_a);
	  tci::copy(ctx,V[site_address_a],A);
	  RankT rank_a = tci::rank(ctx,A);
	  std::vector<BondLabelT> IdxA(rank_a);
	  std::vector<BondLabelT> IdxC(rank_a);
	  for(int k=0; k < bond_idx_a.size(); k++) {
	    auto it_edge_address = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					       bond_idx_a[k]);
	    auto edge_address = std::distance(EdgeIdx.begin(),
					      it_edge_address);
	    std::iota(IdxA.begin(),IdxA.end(),0);
	    std::iota(IdxC.begin(),IdxC.end(),0);
	    IdxA[k] = -1;
	    IdxE[0] = -1;
	    IdxE[1] = k;
	    tci::contract(ctx,A,IdxA,E[edge_address+size_e],IdxE,A,IdxC);
	    ElemT norm = tci::normalize(ctx,A);
	  }
	  tci::copy(ctx,A,AdagA);
	  tci::cplx_conj(ctx,AdagA);
	  std::iota(IdxA.begin(),IdxA.end(),-rank_a);
	  std::iota(IdxC.begin(),IdxC.end(),-rank_a);
	  IdxA[rank_a-1] = 0;
	  IdxC[rank_a-1] = 1;
	  IdxA[target_bond_address_a] = 2;
	  IdxC[target_bond_address_a] = 3;
	  std::vector<BondLabelT> IdxW(4);
	  std::iota(IdxW.begin(),IdxW.end(),0);
	  tci::contract(ctx,A,IdxA,AdagA,IdxC,AdagA,IdxW);
	}

	if( mpi_type == 3 || mpi_type == 1 ) {
	  auto it_site_address_b = std::find(SiteIdx.begin(),SiteIdx.end(),
					     site_b);
	  site_address_b = std::distance(SiteIdx.begin(),it_site_address_b);
	  tci::copy(ctx,V[site_address_b],B);
	  RankT rank_b = tci::rank(ctx,B);
	  std::vector<BondLabelT> IdxB(rank_b);
	  for(int k=0; k < bond_Idx_b.size(); k++) {
	    auto it_edge_address = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					     bond_idx_b[k]);
	    auto edge_address = std::distance(EdgeIdx.begin(),
					      it_edge_address);
	    std::iota(IdxB.begin(),IdxB.end(),0);
	    IdxB[k] = -1;
	    IdxE[0] = -1;
	    IdxE[1] = k;
	    tci::contrat(ctx,B,IdxB,E[edge_address+size_e],IdxE,B);
	    ElemT norm = tci::normalize(ctx,B);
	  }
	  tci::copy(ctx,B,BdagB);
	  tci::cplx_conj(ctx,BdagB);
	  std::vector<BondLabelT> IdxC(rank_b);
	  std::iota(IdxB.begin(),IdxB.end(),-rank_b);
	  std::iota(IdxC.begin(),IdxC.end(),-rank_b);
	  IdxB[rank_b-1] = 0;
	  IdxC[rank_b-1] = 1;
	  IdxB[target_bond_address_b] = 2;
	  IdxC[target_bond_address_c] = 3;
	  std::vector<BondLabelT> IdxW(4);
	  std::iota(IdxW.begin(),IdxW.end(),0);
	  tci::contract(ctx,B,IdxB,BdagB,IdxC,BdagB,IdxW);
	}

	if( mpi_type == 1 ) {
	  MpiSend(ctx,BdagB,mpi_rank_a,comm);
	}

	if( mpi_type == 2 ) {
	  MpiRecv(ctx,BdagB,mpi_rank_b,comm);
	}

	if( mpi_type == 2 || mpi_type == 3 ) {
	  TenT AIA;
	  TenT BIB;
	  TenT S;
	  TenT R;
	  std::vector<std::pair<BondIdxT,BondIdxT>> trace_label(1);
	  trace_label = std::make_pair<BondIdxT,BondIdxT>(0,1);
	  tci::trace(ctx,AdagA,trace_label,AIA);
	  tci::trace(ctx,BdagB,trace_label,BIB);
	  RankT rank_a = tci::rank(ctx,AdagA);
	  std::vector<BondLabelT> IdxA(4);
	  std::vector<BondLabelT> IdxO(4);
	  std::vector<BondLabelT> IdxC(4);
	  IdxA[0] = -1;
	  IdxA[1] = -2;
	  IdxA[2] = 2;
	  IdxA[3] = 3;
	  IdxO[0] = -1;
	  IdxO[1] = 0;
	  IdxO[2] = -2;
	  IdxO[3] = 1;
	  tci::contract(ctx,AdagA,IdxA,O[m],IdxO,AdagA,IdxC);
	  RankT rank_b = tci::rank(ctx,BdagB);
	  std::vector<BondLabelT> IdxB(rank_b);
	  std::iota(IdxA.begin(),IdxA.end(),-rank_a);
	  std::iota(IdxB.begin(),IdxB.end(),-rank_b);
	  std::vector<BondLabelT> IdxR;
	  tci::contract(ctx,AdagA,IdxA,BdagB,IdxB,R,IdxR);
	  std::vector<BondLabelT> IdxAIA(2);
	  std::vector<BondLabelT> IdxBIB(2);
	  std::iota(IdxAIA.begin(),IdxAIA.end(),-2);
	  std::iota(IdxBIB.begin(),IdxBIB.end(),-2);
	  tci::contract(ctx,AIA,IdxAIA,BIB,IdxBIB,S,IdxR);
	  RankT rank_r = tci::rank(ctx,R);
	  CoorsT coors(rank_r,0);
	  ElemT res_r = tci::get_elem(ctx,R,coors);
	  ElemT res_s = tci::get_elem(ctx,W,coors);
	  res[m] = res_r/res_s;
	}
      }
    }
  }
  
}

#endif
