/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/init.h
@brief Initializer for necessary data for tensor product state
*/
#ifndef TNBP_PTNS_FUNCTION_INIT_H
#define TNBP_PTNS_FUNCTION_INIT_H

#include "tnbp/framework/graph.h"

namespace tnbp {  

  /**
     Standard initializer to setup labels and tensors
     for real-space parallelized tensor product state from bond labels
   */
  template <typename TenT>
  void InitTensorProductState(context_handle_t<TenT> & ctx,
		 const std::vector<std::pair<int,int>> & I,
		 const std::vector<int> & PhysicalBondDim,
		 std::vector<TenT> & V,
		 std::vector<int> & SiteIdx,
		 std::map<int,int> & Site_To_MpiRank,
		 std::vector<TenT> & E,
		 std::vector<int> & EdgeIdx,
		 MPI_Comm comm) {
    
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;

    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    auto Site = GetSiteIndexFromBond(I);
    int NumSites = static_cast<int>(Site.size());
    int SiteIdxStart = 0;
    int SiteIdxEnd   = NumSites;
    get_range(mpi_size,mpi_rank,SiteIdxStart,SiteIdxEnd);
    SiteIdx.resize(SiteIdxEnd-SiteIdxStart);
    std::copy(Site.begin() + SiteIdxStart,
	      Site.begin() + SiteIdxEnd,
	      SiteIdx.begin());
    EdgeIdx.resize(0);
    V.resize(SiteIdxEnd-SiteIdxStart);

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
      *itV = tci::assign_from_range<TenT>(ctx,BondDimV,itDataV,
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
    tci::fill(ctx,BondDimE,static_cast<ElemT>(1.0),Eorig);
    E.resize(2*NumEdges);
    for(auto & Em : E) {
      tci::copy(ctx,Eorig,Em);
    }
    std::vector<int> Site_To_MpiRank_Vector(NumSites,0);
    std::vector<int> Site_To_MpiRank_Send(NumSites,0);
    for(const auto & i : SiteIdx) {
      auto it_site_address = std::find(Site.begin(),Site.end(),i);
      auto site_address = std::distance(Site.begin(),it_site_address);
      Site_To_MpiRank_Send[site_address] = mpi_rank;
    }
    MPI_Allreduce(Site_To_MpiRank_Send.data(),
		  Site_To_MpiRank_Vector.data(),
		  NumSites,
		  MPI_INT,
		  MPI_SUM,
		  comm);
    for(size_t i=0; i < NumSites; i++) {
      Site_To_MpiRank[Site[i]] = Site_To_MpiRank_Vector[i];
    }
  }

  /**
     Initializer for TensorProductState class
   */
  template <typename TenT>
  void Init(context_handle_t<TenT> & ctx,
	    const std::vector<int> & pdim,
	    const std::vector<std::pair<int,int>> & I,
	    MPI_Comm comm,
	    TensorProductState<TenT> & W,
	    std::vector<TenT> & E,
	    std::vector<int> & EdgeIdx) {
    W.I_ = I;
    W.comm_ = comm;
    InitTensorProductState(ctx,I,pdim,W.V_,W.SiteIdx_,
			   W.Site_To_MpiRank_,
			   E,EdgeIdx,W.comm_);
  }
}

#endif
