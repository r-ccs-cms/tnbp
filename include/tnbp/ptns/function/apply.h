/// This file is a part of r-ccs-cms/tnbp
/**
   @file tnbp/ptns/function/apply.h
   @brief attach tpo to tps
 */
#ifndef TNBP_PTNS_FUNCTION_APPLY_H
#define TNBP_PTNS_FUNCTION_APPLY_H

namespace tnbp {

  /**
     Attach tensor product operator to tensor product state
   */
  template <typename TenT>
  void AbsorbTPO(context_handle_t<TenT> & ctx,
		 const std::vector<TenT> & O,
		 const std::vector<std::pair<int,int>> & Edge,
		 std::vector<TenT> & V,
		 const std::vector<int> & SiteIdx,
		 const std::map<int,int> & Site_To_MpiRank,
		 std::vector<TenT> & E,
		 const std::vector<int> & EdgeIdx,
		 MPI_Comm comm) {

    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    auto global_sites = GetSiteIndexFromBond(Edge);

    for(int address=0; address < SiteIdx.size(); address++) {
      int site = SiteIdx[address];
      auto it_global_site_address = std::find(global_sites.begin(),
					      global_sites.end(),
					      site);
      auto global_site_address = std::distance(global_sites.begin(),
					       it_global_site_address);
      RankT rank_v = tci::rank(ctx,V[address]);
      RankT rank_o = tci::rank(ctx,O[global_site_address]);
      RankT num_vb = rank_v-1;
      List<BondLabelT> label_v(rank_v);
      List<BondLabelT> label_o(rank_o);
      List<BondLabelT> label_r(rank_v+rank_o-2);
      ShapeT new_shape_v(rank_v);
      ShapeT shape_v = tci::shape(ctx,V[address]);
      ShapeT shape_o = tci::shape(ctx,O[global_site_address]);
      for(RankT k=0; k < num_vb; k++) {
	label_v[k] = static_cast<BondLabelT>(2*k+0);
	label_o[k] = static_cast<BondLabelT>(2*k+1);
	label_r[2*k+0] = static_cast<BondLabelT>(2*k+0);
	label_r[2*k+1] = static_cast<BondLabelT>(2*k+1);
	new_shape_v[k] = static_cast<BondDimT>(shape_v[k]*shape_o[k]);
      }
      label_v[num_vb] = static_cast<BondLabelT>(-1);
      label_o[num_vb] = static_cast<BondLabelT>(-1);
      label_o[num_vb+1] = static_cast<BondLabelT>(2*num_vb+1);
      label_r[2*num_vb] = static_cast<BondLabelT>(2*num_vb+1);
      new_shape_v[num_vb] = shape_o[num_vb];
      tci::contract(ctx,O[global_site_address],label_o,V[address],label_v,
		    V[address],label_r);
      tci::reshape(ctx,V[address],new_shape_v);
    }

    size_t num_e = EdgeIdx.size();
    List<BondLabelT> Idx_I(2);
    List<BondLabelT> Idx_E(2);
    List<BondLabelT> Idx_T(4);
    Idx_E[0] = static_cast<BondLabelT>(0);
    Idx_E[1] = static_cast<BondLabelT>(2);
    Idx_I[0] = static_cast<BondLabelT>(1);
    Idx_I[1] = static_cast<BondLabelT>(3);
    std::iota(Idx_T.begin(),Idx_T.end(),0);

    for(size_t address=0; address < num_e; address++) {
      int site_a = Edge[EdgeIdx[address]].first;
      int site_b = Edge[EdgeIdx[address]].second;
      int mpi_rank_a = Site_To_MpiRank.at(site_a);
      int mpi_rank_b = Site_To_MpiRank.at(site_b);
      
      if( mpi_rank_a == mpi_rank ) {

	auto it_global_site_address_a = std::find(global_sites.begin(),
						  global_sites.end(),
						  site_a);
	auto global_site_address_a = std::distance(global_sites.begin(),
						   it_global_site_address_a);
	std::vector<int> bond = GetSurroundingBondIndex(site_a,Edge);
	
	int target_edge = 0;
	int target_bond = 0;
	for(size_t k=0; k < bond.size(); k++) {
	  if( ( Edge[bond[k]].first == site_a ) &&
	      ( Edge[bond[k]].second == site_b ) ) {
	    target_edge = bond[k];
	    target_bond = k;
	    break;
	  } else if ( (Edge[bond[k]].first == site_b ) &&
		      (Edge[bond[k]].second == site_a ) ) {
	    target_edge = bond[k];
	    target_bond = k;
	    break;
	  }
	}

	auto shape_O = tci::shape(ctx,O[global_site_address_a]);
	auto shape_E = tci::shape(ctx,E[address]);
	ShapeT shape_I(2,shape_O[target_bond]);
	ShapeT shape_N(2,shape_I[0]*shape_E[0]);
	std::vector<ElemT> data_I(shape_O[target_bond]*shape_O[target_bond],
				  static_cast<ElemT>(0.0));
	for(size_t k=0; k < shape_O[target_bond]; k++) {
	  data_I[k+shape_O[target_bond]*k] = static_cast<ElemT>(1.0);
	}
	TenT I;
	auto it_data_I = data_I.begin();
	tci::assign_from_range(ctx,shape_I,it_data_I,
	     [&shape_I](const auto & coor){
	       return coor[0] + shape_I[0] * coor[1];
	     },I);
	TenT T;
	tci::contract(ctx,E[address],Idx_E,I,Idx_I,T,Idx_T);
	tci::reshape(ctx,T,shape_N,E[address]);
	tci::contract(ctx,E[address+num_e],Idx_E,I,Idx_I,T,Idx_T);
	tci::reshape(ctx,T,shape_N,E[address+num_e]);
	
      } else if ( mpi_rank_b == mpi_rank ) {

	auto it_global_site_address_b = std::find(global_sites.begin(),
						  global_sites.end(),
						  site_b);
	auto global_site_address_b = std::distance(global_sites.begin(),
						   it_global_site_address_b);
	
	std::vector<int> bond = GetSurroundingBondIndex(site_b,Edge);

	int target_edge = 0;
	int target_bond = 0;
	for(size_t k=0; k < bond.size(); k++) {
	  if( ( Edge[bond[k]].first == site_b ) &&
	      ( Edge[bond[k]].second == site_a ) ) {
	    target_edge = bond[k];
	    target_bond = k;
	    break;
	  } else if ( (Edge[bond[k]].first == site_a ) &&
		      (Edge[bond[k]].second == site_b ) ) {
	    target_edge = bond[k];
	    target_bond = k;
	    break;
	  }
	}
	
	auto shape_O = tci::shape(ctx,O[global_site_address_b]);
	auto shape_E = tci::shape(ctx,E[address]);
	ShapeT shape_I(2,shape_O[target_bond]);
	ShapeT shape_N(2,shape_I[0]*shape_E[0]);
	std::vector<ElemT> data_I(shape_O[target_bond]*shape_O[target_bond],
				  static_cast<ElemT>(0.0));
	for(size_t k=0; k < shape_O[target_bond]; k++) {
	  data_I[k+shape_O[target_bond]*k] = static_cast<ElemT>(1.0);
	}
	TenT I;
	auto it_data_I = data_I.begin();
	tci::assign_from_range(ctx,shape_I,it_data_I,
	     [&shape_I](const auto & coor){
	       return coor[0] + shape_I[0] * coor[1];
	     },I);
	TenT T;
	tci::contract(ctx,E[address],Idx_E,I,Idx_I,T,Idx_T);
	tci::reshape(ctx,T,shape_N,E[address]);
	tci::contract(ctx,E[address+num_e],Idx_E,I,Idx_I,T,Idx_T);
	tci::reshape(ctx,T,shape_N,E[address+num_e]);
      }
    }
  }
  
  
}

#endif

