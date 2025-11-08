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
			    const std::map<int,int> & Site_To_MpiRank,
			    const std::vector<TenT> & E,
			    const std::vector<int> EdgeIdx,
			    MPI_Comm comm,
			    const std::vector<int> & Site,
			    const std::vector<TenT> & O) {

    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using OrderT = typename tci::tensor_traits<TenT>::order_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t size_e = EdgeIdx.size();

    std::vector<elem_t<TenT>> result(Site.size(),elem_t<TenT>(0.0));

    for(int site_address=0; site_address < Site.size(); site_address++) {
      if( mpi_rank == Site_To_MpiRank.at(Site[site_address]) ) {
	auto it_siteidx_address = std::find(SiteIdx.begin(),
					    SiteIdx.end(),
					    Site[site_address]);
	auto siteidx_address = std::distance(SiteIdx.begin(),
					     it_siteidx_address);
	
	TenT W;
	tci::copy(ctx,V[siteidx_address],W);
	TenT Wdag;
	tci::copy(ctx,V[siteidx_address],Wdag);
	tci::cplx_conj(ctx,Wdag);
	OrderT order_w = tci::order(ctx,W);
	List<BondLabelT> IdxW(order_w);
	List<BondLabelT> IdxC(order_w);
	std::vector<int> BondIdx =
	  GetSurroundingBondIndex(Site[site_address],I);
	for(size_t m=0; m < BondIdx.size(); m++) {
	  auto it_edgeidx_address = std::find(EdgeIdx.begin(),
					      EdgeIdx.end(),
					      BondIdx[m]);
	  auto edgeidx_address = std::distance(EdgeIdx.begin(),
					       it_edgeidx_address);
	  if( I[BondIdx[m]].first == Site[site_address] ) {
	    edgeidx_address += size_e;
	  }
	  List<BondLabelT> IdxF(2);
	  std::iota(IdxW.begin(),IdxW.end(),0);
	  std::iota(IdxC.begin(),IdxC.end(),0);
	  IdxW[m] = static_cast<BondLabelT>(-1);
	  IdxF[0] = static_cast<BondLabelT>(-1);
	  IdxF[1] = static_cast<BondLabelT>(m);
	  tci::contract(ctx,W,IdxW,E[edgeidx_address],IdxF,W,IdxC);
	}

	TenT R;
	List<BondLabelT> IdxO(2);
	std::iota(IdxW.begin(),IdxW.end(),0);
	std::iota(IdxC.begin(),IdxC.end(),0);
	IdxW[order_w-1] = static_cast<BondLabelT>(-1);
	IdxO[0] = static_cast<BondLabelT>(-1);
	IdxO[1] = static_cast<BondLabelT>(order_w-1);
	tci::contract(ctx,W,IdxW,O[site_address],IdxO,R,IdxC);
	std::iota(IdxW.begin(),IdxW.end(),-order_w);
	std::iota(IdxC.begin(),IdxC.end(),-order_w);
	List<BondLabelT> IdxM;
	tci::contract(ctx,W,IdxW,Wdag,IdxC,W,IdxM);
	tci::contract(ctx,R,IdxW,Wdag,IdxC,R,IdxM);
	OrderT order_r = tci::order(ctx,R);
	CoorsT coors(order_r,0);
	ElemT res_r = tci::get_elem(ctx,R,coors);
	ElemT res_w = tci::get_elem(ctx,W,coors);
	result[site_address] = res_r / res_w;
      }
    }

    return result;
  }

  /**
     Function for the measurements of nearest-neighbor two-site local operators
   */
  template <typename TenT>
  std::vector<elem_t<TenT>> Measure(context_handle_t<TenT> & ctx,
			    const std::vector<std::pair<int,int>> & I,
			    const std::vector<TenT> & V,
			    const std::vector<int> & SiteIdx,
			    const std::map<int,int> & Site_To_MpiRank,
			    const std::vector<TenT> & E,
			    const std::vector<int> EdgeIdx,
			    MPI_Comm comm,
			    const std::vector<std::pair<int,int>> & Edge,
			    const std::vector<TenT> & O) {

    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using OrderT = typename tci::tensor_traits<TenT>::order_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t size_e = EdgeIdx.size();

    std::vector<ElemT> result(I.size(),ElemT(0.0));

    for(int m=0; m < Edge.size(); m++) {
      int site_a = Edge[m].first;
      int site_b = Edge[m].second;
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
	auto it_target_bond_address_b = std::find(bond_idx_b.begin(),
						  bond_idx_b.end(),
						  target_edge_address);
	target_bond_address_b = std::distance(bond_idx_b.begin(),
					      it_target_bond_address_b);
	int site_address_a = 0;
	int site_address_b = 0;
	TenT A;
	TenT B;
	TenT AdagA;
	TenT BdagB;
	List<BondLabelT> IdxE(2);

	if( mpi_type == 3 || mpi_type == 2 ) {
	  auto it_site_address_a = std::find(SiteIdx.begin(),
					     SiteIdx.end(),site_a);
	  site_address_a = std::distance(SiteIdx.begin(),it_site_address_a);
	  tci::copy(ctx,V[site_address_a],A);
	  OrderT order_a = tci::order(ctx,A);
	  List<BondLabelT> IdxA(order_a);
	  List<BondLabelT> IdxC(order_a);
	  for(int k=0; k < bond_idx_a.size(); k++) {
	    if( k != target_bond_address_a ) {
	      auto it_edge_address = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					       bond_idx_a[k]);
	      auto edge_address = std::distance(EdgeIdx.begin(),
						it_edge_address);
	      if( I[bond_idx_a[k]].first == site_a ) {
		edge_address += size_e;
	      }
	      std::iota(IdxA.begin(),IdxA.end(),0);
	      std::iota(IdxC.begin(),IdxC.end(),0);
	      IdxA[k] = static_cast<BondLabelT>(-1);
	      IdxE[0] = static_cast<BondLabelT>(-1);
	      IdxE[1] = static_cast<BondLabelT>(k);
	      tci::contract(ctx,A,IdxA,E[edge_address],IdxE,A,IdxC);
	      auto norm = tci::normalize(ctx,A);
	    }
	  }
	  tci::copy(ctx,V[site_address_a],AdagA);
	  tci::cplx_conj(ctx,AdagA);
	  std::iota(IdxA.begin(),IdxA.end(),-order_a);
	  std::iota(IdxC.begin(),IdxC.end(),-order_a);
	  IdxA[order_a-1] = 0;
	  IdxC[order_a-1] = 1;
	  IdxA[target_bond_address_a] = 2;
	  IdxC[target_bond_address_a] = 3;
	  List<BondLabelT> IdxW(4);
	  std::iota(IdxW.begin(),IdxW.end(),0);
	  tci::contract(ctx,A,IdxA,AdagA,IdxC,AdagA,IdxW);
	}

	if( mpi_type == 3 || mpi_type == 1 ) {
	  auto it_site_address_b = std::find(SiteIdx.begin(),SiteIdx.end(),
					     site_b);
	  site_address_b = std::distance(SiteIdx.begin(),it_site_address_b);
	  tci::copy(ctx,V[site_address_b],B);
	  OrderT order_b = tci::order(ctx,B);
	  List<BondLabelT> IdxB(order_b);
	  List<BondLabelT> IdxC(order_b);
	  for(int k=0; k < bond_idx_b.size(); k++) {
	    if( k != target_bond_address_b ) {
	      auto it_edge_address = std::find(EdgeIdx.begin(),
					       EdgeIdx.end(),
					       bond_idx_b[k]);
	      auto edge_address = std::distance(EdgeIdx.begin(),
						it_edge_address);
	      if( I[bond_idx_b[k]].first == site_b ) {
		edge_address += size_e;
	      }
	      std::iota(IdxB.begin(),IdxB.end(),0);
	      std::iota(IdxC.begin(),IdxC.end(),0);
	      IdxB[k] = static_cast<BondLabelT>(-1);
	      IdxE[0] = static_cast<BondLabelT>(-1);
	      IdxE[1] = static_cast<BondLabelT>(k);
	      tci::contract(ctx,B,IdxB,E[edge_address],IdxE,B,IdxC);
	      auto norm = tci::normalize(ctx,B);
	    }
	  }
	  tci::copy(ctx,V[site_address_b],BdagB);
	  tci::cplx_conj(ctx,BdagB);
	  std::iota(IdxB.begin(),IdxB.end(),-order_b);
	  std::iota(IdxC.begin(),IdxC.end(),-order_b);
	  IdxB[order_b-1] = 0;
	  IdxC[order_b-1] = 1;
	  IdxB[target_bond_address_b] = 2;
	  IdxC[target_bond_address_b] = 3;
	  List<BondLabelT> IdxW(4);
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
	  List<Pair<BondIdxT,BondIdxT>> trace_label(1);
	  trace_label[0] = std::make_pair(static_cast<BondIdxT>(0),
				       static_cast<BondIdxT>(1));
	  tci::trace(ctx,AdagA,trace_label,AIA);
	  tci::trace(ctx,BdagB,trace_label,BIB);
	  OrderT order_a = tci::order(ctx,AdagA);
	  List<BondLabelT> IdxA(4);
	  List<BondLabelT> IdxO(4);
	  List<BondLabelT> IdxC(4);
	  IdxA[0] = static_cast<BondLabelT>(-1);
	  IdxA[1] = static_cast<BondLabelT>(-2);
	  IdxA[2] = static_cast<BondLabelT>(2);
	  IdxA[3] = static_cast<BondLabelT>(3);
	  IdxO[0] = static_cast<BondLabelT>(-1);
	  IdxO[1] = static_cast<BondLabelT>(0);
	  IdxO[2] = static_cast<BondLabelT>(-2);
	  IdxO[3] = static_cast<BondLabelT>(1);
	  std::iota(IdxC.begin(),IdxC.end(),0);
	  tci::contract(ctx,AdagA,IdxA,O[m],IdxO,AdagA,IdxC);
	  OrderT order_b = tci::order(ctx,BdagB);
	  List<BondLabelT> IdxB(order_b);
	  std::iota(IdxA.begin(),IdxA.end(),-order_a);
	  std::iota(IdxB.begin(),IdxB.end(),-order_b);
	  List<BondLabelT> IdxR;
	  tci::contract(ctx,AdagA,IdxA,BdagB,IdxB,R,IdxR);
	  List<BondLabelT> IdxAIA(2);
	  List<BondLabelT> IdxBIB(2);
	  std::iota(IdxAIA.begin(),IdxAIA.end(),-2);
	  std::iota(IdxBIB.begin(),IdxBIB.end(),-2);
	  tci::contract(ctx,AIA,IdxAIA,BIB,IdxBIB,S,IdxR);
	  OrderT order_r = tci::order(ctx,R);
	  CoorsT coors(order_r,0);
	  ElemT res_r = tci::get_elem(ctx,R,coors);
	  ElemT res_s = tci::get_elem(ctx,S,coors);
	  result[m] = res_r/res_s;
	}
      }
    }

    return result;
  }
  
}

#endif
