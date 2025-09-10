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
		  real_t<TenT> eps) {
    
  }
  
}

#endif
