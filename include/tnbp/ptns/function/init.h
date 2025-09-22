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
		 std::vector<int> & Site_To_MpiRank,
		 std::vector<TenT> & E,
		 std::vector<int> & EdgeIdx,
		 MPI_Comm comm) {
    
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;

    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

#ifdef INIT_DEBUG
    MPI_Barrier(comm);
    if( mpi_rank == 0 ) {
      std::cout << " tnbp::InitTensorProductState: Start construction of initial TPO " << std::endl;
    }
#endif
    auto Site = GetSiteIndexFromBond(I);
    int NumSites = static_cast<int>(Site.size());
    int SiteIdxStart = 0;
    int SiteIdxEnd   = NumSites;
    get_range(mpi_size,mpi_rank,SiteIdxStart,SiteIdxEnd);

    SiteIdx.resize(SiteIdxEnd-SiteIdxStart);
    std::iota(SiteIdx.begin(),SiteIdx.end(),SiteIdxStart);
    EdgeIdx.resize(0);
    V.resize(SiteIdxEnd-SiteIdxStart);

#ifdef INIT_DEBUG
    MPI_Barrier(comm);
    for(int r=0; r < mpi_size; r++) {
      if( mpi_rank == r ) {
	std::cout << " tnbp::InitTensorProductState: Number of sites handled by "
		  << mpi_rank << "-mpi process = " 
		  << SiteIdxEnd-SiteIdxStart << std::endl;
      }
      MPI_Barrier(comm);
    }
#endif
    
    auto itV = V.begin();
    for(auto const & i : SiteIdx) {
      auto BondIdx = GetSurroundingBondIndex(i,I);
      auto NumBonds = BondIdx.size();
      for(auto const & m : BondIdx) {
	EdgeIdx.push_back(m);
      }
      ShapeT BondDimV(NumBonds+1,1);
      BondDimV[NumBonds] = static_cast<BondDimT>(PhysicalBondDim[i]);
      std::vector<ElemT> DataV(PhysicalBondDim[i],0.0);
      DataV[0] = static_cast<ElemT>(1.0);
      auto itDataV = DataV.begin();
      tci::assign_from_container(ctx,BondDimV,itDataV,
				 [](const CoorsT & coors) {
				   return coors[0];},*itV);
      itV++;
    }
    std::sort(EdgeIdx.begin(),EdgeIdx.end());
    EdgeIdx.erase(std::unique(EdgeIdx.begin(),EdgeIdx.end()),
		  EdgeIdx.end());
#ifdef INIT_DEBUG
    MPI_Barrier(comm);
    for(int r=0; r < mpi_size; r++) {
      if( mpi_rank == r ) {
	std::cout << " tnbp::InitTensorProductState: Number of edges handled by "
		  << mpi_rank << "-mpi process = "
		  << EdgeIdx.size() << std::endl;
      }
      MPI_Barrier(comm);
    }
#endif
    
    int NumEdges = static_cast<int>(EdgeIdx.size());
    TenT Eorig;
    ShapeT BondDimE(2,1);
    tci::fill(ctx,BondDimE,static_cast<ElemT>(1.0),Eorig);
    E.resize(2*NumEdges);
    for(auto & Em : E) {
      tci::copy(ctx,Eorig,Em);
    }
    Site_To_MpiRank.resize(NumSites,0);
    std::vector<int> Site_To_MpiRank_Send(NumSites,0);
    for(const auto & i : SiteIdx) {
      Site_To_MpiRank_Send[i] = mpi_rank;
    }
    MPI_Allreduce(Site_To_MpiRank_Send.data(),
		  Site_To_MpiRank.data(),
		  NumSites,
		  MPI_INT,
		  MPI_SUM,
		  comm);
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
