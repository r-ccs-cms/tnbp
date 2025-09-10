/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/ptns/function/init.h
@brief Initializer for necessary data for tensor product state
*/
#ifndef TNBP_PTNS_FUNCTION_INIT_H
#define TNBP_PTNS_FUNCTION_INIT_H

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

    int mpi_size; MPI_Comm_size(comm,&mpi_size);
    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

    auto Site = GetSiteIndexFromBond(I);
    int NumSites = static_cast<int>(Site.size());
    int SiteIdxStart = 0;
    int SiteIdxEnd   = NumSites;
    get_range(mpi_size,mpi_rank,SiteIdxStart,SiteIdxEnd);

    SiteIdx.resize(SiteIdxEnd-SiteIdxStart);
    std::iota(SiteIdx.begin(),SiteIdx.end(),SiteIdxStart);
    EdgeIdx.resize(0);
    V.resize(SiteIdxEnd-SiteIdxStart);
    auto itV = V.begin();
    for(auto const & i : SiteIdx) {
      auto BondIdx = GetSurroundingBondIndex(i,I);
      auto NumBonds = BondIdx.size();
      for(auto const & m : BondIdx) {
	EdgeIdx.push_back(m);
      }
      std::vector<int> BondDimV(NumBonds+1,1);
      BondDimV[NumBonds] = PhysicalBondDim[i];
      std::vector<ElemT> DataV(PhysicalBondDim[i],0.0);
      DataV[0] = ElemT(1.0);
      tci::allocate(ctx,BondDimV,*itV);
      auto itDataV = DataV.begin();
      tci::for_each(ctx,*itV,[&itDataV](ElemT & elem) {
	elem = *itDataV++;
      });
      itV++;
    }
    std::sort(EdgeIdx.begin(),EdgeIdx.end());
    EdgeIdx.erase(std::unique(EdgeIdx.begin(),EdgeIdx.end()),EdgeIdx.end());
    int NumRankEdges = static_cast<int>(EdgeIdx.size());
    TenT Eorg;
    tci::allocate(ctx,{1,1},Eorg);
    tci::for_each(ctx,Eorg,[](ElemT & elem) { elem = static_cast<ElemT>(1.0); });
    E.resize(2*NumRankEdges,Eorg);
    Site_To_MpiRank(NumSites,0);
    std::vector<int> Site_To_MpiRank_Send(NumSites,0);
    for(auto const & i : SiteIdx) {
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
  friend void Init(context_handle_t<TenT> & ctx,
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
			   E,EdgeIdx,W.comm_)
  }
}

#endif
