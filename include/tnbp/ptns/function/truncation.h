/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/truncation.h
@brief truncation based on the messenger tensors
*/
#ifndef TNBP_PTNS_FUNCTION_TRUNCATION_H
#define TNBP_PTNS_FUNCTION_TRUNCATION_H

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

    int mpi_rank; MPI_Comm_rank(comm,mpi_rank);
    int mpi_size; MPI_Comm_size(comm,mpi_size);
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
	    target_edge == bond_a[k];
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
	  SquareRootAndInverse(E[edge_address],Ra,Sa,eps);
	  SquareRootAndInverse(E[edge_address+num_e],Rb,Sb,eps);
	} else {
	  SquareRootAndInverse(E[edge_address+num_e],Ra,Sa,eps);
	  SquareRootAndInverse(E[edge_address],Rb,Sb,eps);
	}

	TenT T;
	std::vector<BondLabelT> IdxRa(2);
	std::vector<BondLabelT> IdxRb(2);
	IdxRa[0] = -1;
	IdxRa[1] = 0;
	IdxRb[0] = -1;
	IdxRb[1] = 1;
	tci::contract(ctx,Ra,IdxRa,Rb,IdxRb,T);
	auto norm_t = tci::normalize(ctx,T);

	TenT X;
	TenT Y;
	RealTenT Zdiag;
	RankT num_rows = 1;
	RealT trunc_err;
	BondDimT chi_min = 1;
	BondDimT chi_max = max_dim;
	tci::trunc_svd(ctx,T,num_rows,X,Zdiag,Y,
		       trunc_err,chi_min,chi_max,err,eps);
	RealTenT Pdiag = tci::copy(ctx,Zdiag);
	tci::for_each(ctx,Pdiag,[](RealT & elem) {
	  elem = std::sqrt(elem); });
	ShapeT shape_z = tci::shape(ctx,Zdiag);
	std::vector<ElemT> data_z(shape_z[0]);
	auto it_data_z = data_z.begin();
	tci::for_each(ctx,Zdiag,[&it_data_z](const RealT & elem) {
	  *it_data_z++ = static_cast<ElemT>(elem); });
	std::vector<ElemT> data_p(shape_z[0]);
	auto it_data_p = data_p.begin();
	tci::for_each(ctx,Pdiag,[&it_data_p](const RealT & elem) {
	  *it_data_p++ = static_cast<ElemT>(elem)});
	
	TenT Z; tci::allocate(ctx,shape_z,Z);
	TenT P; tci::allocate(ctx,shape_z,P);
	it_data_z = data_z.begin();
	it_data_p = data_p.begin();
	tci::for_each(ctx,Z,[&it_data_z](ElemT & elem) {
	  elem = *it_data_z++; });
	tci::for_each(ctx,P,[&it_data_p](ElemT & elem) {
	  elem = *it_data_p++; });
	
	tci::diag(ctx,Z);
	tci::diag(ctx,P);
	tci::copy(ctx,Z,E[edge_address]);

	std::vector<BondLabelT> IdxU(2);
	std::vector<BondLabelT> IdxS(2);
	std::vector<BondLabelT> IdxP(2);
	std::vector<BondLabelT> IdxT(2);

	if( mpi_type == 2 || mpi_type == 3 ) {
	  IdxS[0] = -1;
	  IdxS[1] = 0;
	  IdxU[0] = -1;
	  IdxU[1] = 1;
	  tci::contract(Sa,IdxS,X,IdxU,T);
	  IdxT[0] = 0;
	  IdxT[1] = -1;
	  IdxP[0] = -1;
	  IdxP[1] = 1;
	  tci::contract(T,IdxT,P,IdxP,T);
	  auto it_site_address_a = std::find(SiteIdx.begin(),SiteIdx.end(),
					site_a);
	  site_address_a = std::distance(SiteIdx.begin(),it_site_address_a);

	  RankT rank_a = tci::rank(ctx,V[site_address_a]);
	  std::vector<BondLabelT> IdxA(rank_a);
	  std::iota(IdxA.begin(),IdxA.end(),0);
	  IdxA[bond_address_a] = -1;
	  IdxT[0] = -1;
	  IdxT[1] = bond_address_a;
	  tci::contract(V[site_address_a],IdxA,T,IdxT,
			V[site_address_a]);
	  auto norm_a = tci::normalize(ctx,V[site_address_a]);
	}

	if( mpi_type == 1 || mpi_type == 3 ) {
	  IdxU[0] = 0;
	  IdxU[1] = -1;
	  IdxS[0] = -1;
	  IdxS[1] = 1;
	  tci::contract(Y,IdxU,Sb,IdxS,T);
	  IdxP[0] = 0;
	  IdxP[1] = -1;
	  IdxT[0] = -1;
	  IdxT[1] = 1;
	  tci::contract(P,IdxP,T,IdxT,T);
	  auto it_site_address_b = std::find(SiteIdx.begin(),SiteIdx.end(),
					     site_b);
	  site_address_b = std::distance(SiteIdx.begin(),it_site_address_b);
	  RankT rank_b = tci::rank(ctx,V[site_address_b]);
	  std::vector<BondLabelT> IdxB(rank_b);
	  std::iota(IdxB.begin(),IdxB.end(),0);
	  IdxB[bond_address_b] = -1;
	  IdxT[0] = bond_address_b;
	  IdxT[1] = -1;
	  tci::contract(T,IdxT,V[site_address_b],IdxB,
			V[site_address_b]);
	  auto norm_b = tci::normalize(ctx,V[site_address_b]);
	}
      }
    }

    
  }
  
}

#endif
