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
    TensorProductState(const std::vector<TenT> & V,
		       const std::vector<std::pair<int,int>> & I,
		       const std::vector<int> SiteIdx,
		       const std::vector<int> Site_To_MpiRank,
		       MPI_Comm comm) :
      V_(V), I_(I), SiteIdx_(SiteIdx),
      Site_To_MpiRank_(Site_To_MpiRank), comm_(comm) {
      MPI_Comm_size(comm_,&mpi_size_);
      MPI_Comm_rank(comm_,&mpi_rank_);
    }

    /**
       Constructor from other TensorProductState
     */
    TensorProductState(const TensorProductState & other) :
      V_(other.V_), I_(other.I_), SiteIdx_(other.SiteIdx_),
      Site_To_MpiRank_(other.Site_To_MpiRank_),
      comm_(other.comm_), mpi_size_(other.mpi_size_),
      mpi_rank_(other.mpi_rank_) {}

    /**
       Deconstructor
     */
    ~TensorProductState() {}

    /**
       Copy operator
     */
    TensorProductState & operator = (const TensorProductState & other) {
      if( this != &other ) {
	copy(other);
      }
      return *this;
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
    std::vector<int> Site_To_MpiRank_;

    MPI_Comm comm_;
    int mpi_master_;
    int mpi_size_;
    int mpi_rank_;

    void copy(const TensorProductState & other) {
      this->V_ = other.V_;
      this->I_ = other.I_;
      this->SiteIdx_ = other.SiteIdx_;
      this->Site_To_MpiRank_ = other.Site_To_MpiRank_;
      this->comm_ = other.comm_;
      this->mpi_size_ = other.mpi_size_;
      this->mpi_rank_ = other.mpi_rank_;
    }
  };
  
}

#endif
