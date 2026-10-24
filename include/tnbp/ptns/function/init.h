/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/init.h
@brief Initializer for necessary data for tensor product state
*/
#ifndef TNBP_PTNS_FUNCTION_INIT_H
#define TNBP_PTNS_FUNCTION_INIT_H

#include "tnbp/framework/mpiutility.h"
#include "tnbp/framework/typedef.h"
#include "tnbp/ptns/tps.h"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tnbp/framework/graph.h"
#include "tnbp/ptns/function/message_gauge.h"

namespace tnbp {  

  /**
     Standard initializer to setup labels and tensors
     for real-space parallelized tensor product state from bond labels.
     The caller supplies the complete owner map; this function never partitions sites.
   */
  template <typename TenT>
  void InitTensorProductState(context_handle_t<TenT> & ctx,
		 const std::vector<std::pair<int,int>> & I,
		 const std::vector<int> & PhysicalBondDim,
		 std::vector<TenT> & V,
		 std::vector<int> & SiteIdx,
		 const std::map<int,int> & Site_To_MpiRank,
		 std::vector<TenT> & E,
		 std::vector<int> & EdgeIdx,
		 MPI_Comm comm) {
    
    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;
    using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
    using BondDimT = typename tcapi::tensor_traits<TenT>::bond_dim_t;

    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    auto Site = GetSiteIndexFromBond(I);
    ValidateSiteOwners(Site,Site_To_MpiRank,comm);
    SiteIdx.clear();
    for(auto site:Site) if(Site_To_MpiRank.at(site)==mpi_rank) SiteIdx.push_back(site);
    EdgeIdx.clear();
    V.resize(SiteIdx.size());

    auto itV = V.begin();
    for(auto const & i : SiteIdx) {
      auto BondIdx = GetSurroundingBondIndex(i,I);
      auto NumBonds = BondIdx.size();
      auto it_site_address = std::find(Site.begin(),Site.end(),i);
      auto site_address = std::distance(Site.begin(),it_site_address);
      for(auto const & m : BondIdx) {
	EdgeIdx.push_back(m);
      }
      ShapeT BondDimV(NumBonds+1,1);
      BondDimV[NumBonds] = static_cast<BondDimT>(PhysicalBondDim[site_address]);
      std::vector<ElemT> DataV(PhysicalBondDim[site_address],0.0);
      DataV[0] = static_cast<ElemT>(1.0);
      auto itDataV = DataV.begin();
      *itV = tcapi::assign_from_range<TenT>(ctx,BondDimV,itDataV,
		    [&BondDimV](const CoorsT & coors) {
		      return address_from_coor(BondDimV,coors); });
      itV++;
    }
    std::sort(EdgeIdx.begin(),EdgeIdx.end());
    EdgeIdx.erase(std::unique(EdgeIdx.begin(),EdgeIdx.end()),
		  EdgeIdx.end());
    
    int NumEdges = static_cast<int>(EdgeIdx.size());
    TenT Eorig;
    ShapeT BondDimE(2,1);
    Eorig = tcapi::fill<TenT>(ctx,BondDimE,static_cast<ElemT>(1.0));
    NormalizeMessage(ctx,Eorig,comm);
    E.resize(2*NumEdges);
    for(auto & Em : E) {
      Em = tcapi::copy(ctx,Eorig);
    }

  }

  /**
     Initializer for TensorProductState class
   */
  template <typename TenT>
  void Init(context_handle_t<TenT> & ctx,
	    const std::vector<int> & pdim,
	    const std::vector<std::pair<int,int>> & I,
            const std::map<int,int> & owners,
	    MPI_Comm comm,
	    TensorProductState<TenT> & W,
	    std::vector<TenT> & E,
	    std::vector<int> & EdgeIdx) {
    W.Site_To_MpiRank_ = owners;
    W.I_ = I;
    W.comm_ = comm;
    InitTensorProductState(ctx,I,pdim,W.V_,W.SiteIdx_,
			   W.Site_To_MpiRank_,
			   E,EdgeIdx,W.comm_);
  }

  /**
     @brief Initializer for the messenger tensor for edges
   */
  template <typename TenT, typename IntT>
  void init_edge_messenger_tensors(context_handle_t<TenT> & ctx,
				   const std::vector<std::pair<IntT,IntT>> & edges,
				   const std::vector<TenT> & V,
				   const std::vector<IntT> & site_idx,
				   const std::map<IntT,IntT> & site_to_mpi_rank,
				   std::vector<TenT> & E,
				   std::vector<IntT> & edge_idx,
				   MPI_Comm comm) {
    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    
    edge_idx.resize(0);
    E.resize(0);
    for(auto const & site : site_idx) {
      auto bonds = GetSurroundingBondIndex(site,edges);
      for(auto const & m : bonds) {
	edge_idx.push_back(m);
      }
    }
    std::sort(edge_idx.begin(),edge_idx.end());
    edge_idx.erase(std::unique(edge_idx.begin(),edge_idx.end()),edge_idx.end());
    auto num_edges = edge_idx.size();
    E.resize(2*num_edges);

    // define bond dim from the site_a vertex tensor
    for(std::size_t edge_adrs=0;
	edge_adrs < num_edges;
	edge_adrs++) {
      auto site_a = edges[edge_idx[edge_adrs]].first;
      auto site_b = edges[edge_idx[edge_adrs]].second;
      int mpi_type = 0;
      auto mpi_rank_a = site_to_mpi_rank.at(site_a);
      auto mpi_rank_b = site_to_mpi_rank.at(site_b);
      if( mpi_rank_a == mpi_rank ) {
	mpi_type += 2;
      }
      if( mpi_rank_b == mpi_rank ) {
	mpi_type += 1;
      }
      if( mpi_type == 2 || mpi_type == 3 ) {
	auto it_site_adrs = std::find(site_idx.begin(),site_idx.end(),
				      site_a);
	auto site_adrs = std::distance(site_idx.begin(),it_site_adrs);
	auto shape_v = tcapi::shape(ctx,V[site_adrs]);
	auto bond = GetSurroundingBondIndex(site_a,edges);
	auto it_target_bond = std::find(bond.begin(),bond.end(),edge_idx[edge_adrs]);
	auto target_bond = std::distance(bond.begin(),it_target_bond);
	E[edge_adrs] = tcapi::eye<TenT>(ctx,shape_v[target_bond]);
	E[edge_adrs+num_edges] = tcapi::eye<TenT>(ctx,shape_v[target_bond]);
      } else if ( mpi_type == 1 ) {
	auto it_site_adrs = std::find(site_idx.begin(),site_idx.end(),
				      site_b);
	auto site_adrs = std::distance(site_idx.begin(),it_site_adrs);
	auto shape_v = tcapi::shape(ctx,V[site_adrs]);
	auto bond = GetSurroundingBondIndex(site_b,edges);
	auto it_target_bond = std::find(bond.begin(),bond.end(),edge_idx[edge_adrs]);
	auto target_bond = std::distance(bond.begin(),it_target_bond);
	E[edge_adrs] = tcapi::eye<TenT>(ctx,shape_v[target_bond]);
	E[edge_adrs+num_edges] = tcapi::eye<TenT>(ctx,shape_v[target_bond]);
      }
    }
    NormalizeMessages(ctx,E,comm);
  }
  
}

#endif
