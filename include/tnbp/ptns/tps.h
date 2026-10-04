// This file is a part of r-ccs-cms/tnbp
/**
@file ptns/tps.h
@brief The header file to define the tensor product state
*/
#ifndef TNBP_PTNS_TPS_H
#define TNBP_PTNS_TPS_H

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/mpiutility.h"

namespace tnbp {

  template <typename TenT>
  class TensorProductState {

  public:

    /**
       Default constructor of TensorProductState
     */
    TensorProductState() : V_(), I_(), SiteIdx_(), Site_To_MpiRank_() {}

    /**
       Constructor with specifying all data
     */
    TensorProductState(context_handle_t<TenT> & ctx,
                       const std::vector<TenT> & V,
		       const std::vector<std::pair<int,int>> & I,
		       const std::vector<int> SiteIdx,
		       const std::map<int,int> Site_To_MpiRank,
		       MPI_Comm comm) :
      I_(I), SiteIdx_(SiteIdx),
      Site_To_MpiRank_(Site_To_MpiRank), comm_(comm) {
      V_.reserve(V.size());
      for (const auto& tensor : V) V_.push_back(tcapi::copy(ctx,tensor));
      MPI_Comm_size(comm_,&mpi_size_);
      MPI_Comm_rank(comm_,&mpi_rank_);
    }

    TensorProductState(const TensorProductState&) = delete;
    TensorProductState& operator=(const TensorProductState&) = delete;
    TensorProductState(TensorProductState&&) = default;
    TensorProductState& operator=(TensorProductState&&) = default;

    TensorProductState copy(context_handle_t<TenT>& ctx) const {
      TensorProductState result;
      result.V_.reserve(V_.size());
      for (const auto& tensor : V_) result.V_.push_back(tcapi::copy(ctx,tensor));
      result.I_ = I_;
      result.SiteIdx_ = SiteIdx_;
      result.Site_To_MpiRank_ = Site_To_MpiRank_;
      result.comm_ = comm_;
      result.mpi_master_ = mpi_master_;
      result.mpi_size_ = mpi_size_;
      result.mpi_rank_ = mpi_rank_;
      return result;
    }

    /**
       Number of vertex tensors
     */
    size_t NumV() { return V_.size(); }

    /**
       Site Index
     */
    int Q(size_t i) { return SiteIdx_[i]; }

    /**
       Initializer function
     */
    template <typename TenT_>
    friend void Init(context_handle_t<TenT_> & ctx,
		     const std::vector<int> & pdim,
		     const std::vector<std::pair<int,int>> & I,
		     MPI_Comm comm,
		     TensorProductState<TenT_> & W,
		     std::vector<TenT_> & E,
		     std::vector<int> & EdgeIdx);

    /**
       Function to calculate messenger tensors for belief propagation
     */
    template <typename TenT_>
    friend void BeliefPropagation(context_handle_t<TenT_> & ctx,
			     const TensorProductState<TenT_> & W,
			     const std::vector<std::pair<int,int>> & J,
			     const std::vector<TenT_> & E,
			     const std::vector<int> & EdgeIdx,
			     std::vector<TenT_> & F);

    /**
       Function to evaluate the error from bp condition
     */
    template <typename TenT_>
    friend void BeliefPropagationCondition(context_handle_t<TenT_> & ctx,
				 const TensorProductState<TenT_> & W,
				 const std::vector<TenT_> & E,
				 const std::vector<int> & EdgeIdx,
				 real_t<TenT_> & result);

  private:

    std::vector<TenT> V_;
    std::vector<std::pair<int,int>> I_;
    std::vector<int> SiteIdx_;
    std::map<int,int> Site_To_MpiRank_;

    MPI_Comm comm_ = MPI_COMM_NULL;
    int mpi_master_ = 0;
    int mpi_size_ = 0;
    int mpi_rank_ = 0;


  };
  
}

#endif
