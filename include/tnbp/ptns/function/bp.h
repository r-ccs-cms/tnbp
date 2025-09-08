/// This file is a part of r-ccs-cms/tnbp
/**
   @file tnbp/ptns/function/bp.h
   @brief Belief propagation methods
*/
#ifndef TNBP_PTNS_FUNCTION_BP_H
#define TNBP_PTNS_FUNCTION_BP_H

namespace tnbp {

/**
Belief propagation step for self-consistent iteration
*/
  template <typename TenT>
  void BeliefPropagation(context_handle_t<TenT> & ctx,
			 const std::vector<std::pair<int,int>> & I,
			 const std::vector<TenT> & V,
			 const std::vector<int> & SiteIdx,
			 const std::vector<int> & Site_To_MpiRank,
			 const std::vector<std::pair<int,int>> & J,
			 const std::vector<TenT> & E,
			 const std::vector<int> & EdgeIdx,
			 MPI_Comm comm,
			 std::vector<TenT> & F) {

    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using BondLabelT  = typename tci::tensor_traits<Tent>::bond_label_t;

    int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
    int mpi_size; MPI_Comm_size(comm,&mpi_size);

    size_t size_e = EdgeIdx.size();
    F.resize(2*size_e);

    for(int m=0; m < J.size(); m++) {
      int SiteA = J[m].first;
      int SiteB = J[m].second;
      int MpiRankA = Site_To_MpiRank[SiteA];
      int MpiRankB = Site_To_MpiRank[SiteB];
      int type = 0;

      if( MpiRankA == mpi_rank ) {
	type += 2;
      }
      if( MpiRankB == mpi_rank ) {
	type += 1;
      }

      if( type > 0 ) {
	std::vector<int> EdgeIdxA = GetSurroundingBondIndex(SiteA,I);
	std::vector<int> EdgeIdxB = GetSurroundingBondIndex(SiteB,I);
	int TargetEdgeIdx = 0;
	int BondIdxA = 0;
	int BondIdxB = 0;

	for(int k=0; k < EdgeA.size(); k++) {
	  if( (I[EdgeIdxA[k]].first == SiteA)
	      && ( I[EdgeIdxA[k]].second == SiteB ) ) {
	    TargetEdgeIdx = EdgeIdxA[k];
	    BondIdxA = k;
	  } else if ( (I[EdgeIdxA[k]].first == SiteB)
		      && ( I[EdgeIdxA[k]].second == SiteA) ) {
	    TargetEdgeIdx = EdgeIdxA[k];
	    BondIdxA = k;
	  }
	}

	auto itBondIdxB = std::find(EdgeIdxB.begin(),EdgeIdxB.end(),
				    TargetEdgeIdx);
	int BondIdxB = std::distance(EdgeIdxB.begin(),itBondIdxB);

	int AddressA;
	int AddressB;
	TenT AdagA;
	TenT BdagB;

	if( type == 3 || type == 2 ) {
	  auto itAddressA = std::find(SiteIdx.begin(),SiteIdx.end(),SiteA);
	  AddressA = std::distance(SiteIdx.begin(),itAddressA);
	  TenT T = V[AddressA];
	  auto RankA = tci::rank(ctx,T);
	  std::vector<BondLabelT> IdxA(RankA);
	  std::vector<BondLabelT> IdxE(2);

	  for(int k=0; k < EdgeIdxA.size(); k++) {
	    auto itMpA = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				   EdgeIdxA[k]);
	    int MpA = std::distance(EdgeIdx.begin(),itMpA);
	    if( k != BondIdxA ) {
	      std::iota(IdxA.begin(),IdxA.end(),0);
	      IdxA[k] = -1;
	      IdxE[0] = -1;
	      IdxE[1] = k;
	      if( I[EdgeIdxA[k]].first == SiteA ) {
		tci::contract(ctx,T,IdxA,E[MpA+size_e],IdxE,T);
	      } else {
		tci::contract(ctx,T,IdxA,E[MpA],IdxE,T);
	      }
	    }
	    ElemT norm_a = tci::normalize(ctx,T);
	  }
	  AdagA = V[AddressA];
	  tci::cplx_conj(ctx,AdagA);
	  
	  std::vector<BondLabelT> IdxT(RankA);
	  std::iota(IdxA.begin(),IdxA.end(),-RankA);
	  std::iota(IdxT.begin(),IdxT.end(),-RankA);
	  IdxT[BondIdxA] = 0;
	  IdxA[BondIdxA] = 1;
	  tci::contract(ctx,T,IdxT,AdagA,IdxA,AdagA);
	  ElemT norm_aa = tci::normalize(ctx,AdagA);
	  auto itMpT = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				 TargetEdgeIdx);
	  int MpT = std::distance(EdgeIdx.begin(),itMpT);
	  if( I[TargetEdgeIdx].first == SiteA ) {
	    F[MpT] = AdagA;
	  } else {
	    F[MpT+size_e] = AdagA;
	  }
	  
	}

	if( type == 3 || type == 1 ) {

	  auto itAddressB = std::find(SiteIdx.begin(),SiteIdx.end(),SiteB);
	  AddressB = std::distance(SiteIdx.begin(),itAddressB);
	  TenT T = V[AddressB];
	  auto RankB = tci::rank(ctx,T);
	  std::vector<BondLabelT> IdxB(RankB);
	  std::vector<BondLabelT> IdxE(2);

	  for(int k=0; k < EdgeIdxB.size(); k++) {
	    auto itMpA = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				   EdgeIdxB[k]);
	    int MpB = std::distance(EdgeIdx.begin(),itMpB);
	    if( k != BondIdxB ) {
	      std::iota(IdxB.begin(),IdxB.end(),0);
	      IdxB[k] = -1;
	      IdxE[0] = -1;
	      IdxE[1] = k;
	      if( I[EdgeIdxB[k]].first == SiteB ) {
		tci::contract(ctx,T,IdxB,E[MpB+size_e],IdxE,T);
	      } else {
		tci::contract(ctx,T,IdxB,E[MpB],IdxE,T);
	      }
	    }
	    ElemT norm_b = tci::normalize(ctx,T);
	  }
	  BdagB = V[AddressB];
	  tci::cplx_conj(ctx,BdagB);
	  
	  std::vector<BondLabelT> IdxT(RankB);
	  std::iota(IdxB.begin(),IdxB.end(),-RankB);
	  std::iota(IdxT.begin(),IdxT.end(),-RankB);
	  IdxT[BondIdxB] = 0;
	  IdxA[BondIdxB] = 1;
	  tci::contract(ctx,T,IdxT,BdagB,IdxB,BdagB);
	  ElemT norm_bb = tci::normalize(ctx,BdagB);
	  auto itMpT = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				 TargetEdgeIdx);
	  int MpT = std::distance(EdgeIdx.begin(),itMpT);
	  if( I[TargetEdgeIdx].first == SiteB ) {
	    F[MpT] = BdagB;
	  } else {
	    F[MpT+size_e] = BdagB;
	  }
	  
	}

	if( type == 1 ) {
	  MpiSend(&BdagB,MpiRankA,comm);
	}
	if( type == 2 ) {
	  MpiRedv(&BdagB,MpiRankB,comm);
	  auto itMpT = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				 TargetEdgeIdx);
	  int MpT = std::distance(EdgeIdx.begin(),itMpT);
	  if( I[TargetEdgeIdx].first == SiteB ) {
	    F[MpT] = BdagB;
	  } else {
	    F[MpT+size_e] = BdagB;
	  }
	}

	if( type == 2 ) {
	  MpiSend(AdagA,MpiRankB,comm);
	}
	if( type == 1 ) {
	  MpiRecv(AdagA,MpiRankA,comm);
	  auto itMpT = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				 TargetEdgeIdx);
	  int MpT = std::distance(EdgeIdx.begin(),itMpT);
	  if( I[TargetEdgeIdx].first == SiteA ) {
	    F[MpT] = AdagA;
	  } else {
	    F[MpT+size_e] = AdagA;
	  }
	}
	
      } // end if( type > 0 )
    } // end for(int m=0; m < J.size(); m++)
  }

  /**
     Function to check belief propagation condition
   */
  template <typename TenT, typename RealT>
  void BeliefPropagationCondition(context_handle_t<TenT> & ctx,
				  const std::vector<std::pair<int,int>> & I,
				  const std::vector<TenT> & V,
				  const std::vector<int> & SiteIdx,
				  const std::vector<int> & Site_To_MpiRank,
				  const std::vector<TenT> & E,
				  const std::vector<int> & EdgeIdx,
				  MPI_Comm comm,
				  RealT & result) {
    
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using BondLabelT  = typename tci::tensor_traits<TenT>::bond_label_t;

    size_t num_v = SiteIdx.size();
    size_t num_e = EdgeIdx.size();
    size_t num_total_edges = I.size();
    RealT result_volume = 1.0/(2.0*total_edges);
    RealT result_rank = 0.0;
    result = 0.0;

    for(size_t address=0; address < SiteIdx.size(); address++) {
      std::vector<int> BondIdx = GetSurroundingBondIdx(SiteIdx[address],I);
      for(size_t k=0; k < BondIdx.size(); k++) {
	TenT T = V[address];
	auto RankV = tci::rank(ctx,T);
	for(size_t l=0; l < BondIdx.size(); l++) {
	  if( l != k ) {
	    auto itEdgeAddress = std::find(EdgeIdx.begin(),EdgeIdx.end(),
					   BondIdx[l]);
	    int EdgeAddress = std::distance(EdgeIdx.begin(),itEdgeAddress);
	    if( I[BondIdx[l]].first == SiteIdx[address] ) {
	      EdgeAddress += num_e;
	    }
	    TenT F = E[EdgeAddress];
	    std::vector<BondLabelT> IdxE(2);
	    std::vector<BondLabelT> IdxV(RankV);
	    std::iota(IdxV.begin(),IdxV.end(),0);
	    IdxV[l] = -1;
	    IdxE[0] = -1;
	    IdxE[1] = l;
	    tci::contract(ctx,T,IdxV,F,IdxE,T);
	  }
	} // end for(size_t l=0; l < BondIdx.size(); l++)
	
	TenT C = V[address];
	tci::cplx_conj(ctx,C);
	std::vector<BondLabelT> IdxV(RankV);
	std::vector<BondLabelT> IdxC(RankV);
	std::iota(IdxV.begin(),IdxV.end(),-RankV);
	std::iota(IdxC.begin(),IdxC.end(),-RankV);
	IdxV[k] = 0;
	IdxC[k] = 1;
	tci::contract(ctx,T,IdxV,C,IdxC,T);
	ElemT norm_v = tci::normalize(ctx,T);

	auto itEdgeAddress = std::find(EdgeIdx.begin(),EdgeIdx.end(),
				       BondIdx[k]);
	int EdgeAddress = std::distance(EdgeIdx.begin(),itEdgeAddress);
	if( I[BondIdx[k]].second == SiteIdx[address] ) {
	  EdgeAddress += num_e;
	}
	TenT F = E[EdgeAddress];
	ElemT norm_f = tci::normalize(ctx,F);

	TenT D;
	tci::linear_combine(ctx,{T,F},{ElemT(1.0),ElemT(-1.0)},D);
	RealT NormD = tci::norm(ctx,D);
	result_rank += NormD * result_volume;
	
      } // end for(size_t k=0; k < BondIdx.size(); k++)
    } // end for(int address=0; address < SiteIdx.size(); address++)
    MPI_Allreduce(&result_rank,&result,1,MPI_DOUBLE,MPI_SUM,comm);
  }

  

}

#endif

