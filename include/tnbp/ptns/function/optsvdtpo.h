/// This file is a part of r-ccs-cms/tnbp
/**
   @file tnbp/ptns/function/optsvdtpo.h
   @brief Simple function to optimize tpo by svd
 */
#ifndef TNBP_PTNS_FUNCTION_OPTSVDTPO_H
#define TNBP_PTNS_FUNCTION_OPTSVDTPO_H

#include "tnbp/framework/typedef.h"
#include "tnbp/framework/graph.h"
#include "tnbp/framework/mpiutility.h"

namespace tnbp {

  /**
   * @brief Optimizes the tensor product operator (TPO) @p O by contracting and refactorizing tensors via SVD, truncating small singular values below a given threshold.
   *
   * @details
   * Under the TCI context @p ctx, this function processes each tensor pair specified in @p edges.
   * For each bond pair (i, j):
   * 1) The tensors associated with the bond (typically O[i] and O[j]) are temporarily contracted over the shared indices.
   * 2) The resulting tensor is reshaped into a two-dimensional form suitable for singular value decomposition (SVD).
   * 3) Singular values smaller than @p eps are discarded (hard cutoff), and the truncated SVD factors are reshaped and reassigned to the corresponding tensors.
   *
   * This effectively reduces the bond dimension of the operator while maintaining the essential correlations between tensors.
   * It is mainly used for compressing Matrix Product Operators (MPOs) or Tensor Product Operators (TPOs) prior to simulation or optimization.
   * The truncation threshold @p eps controls the trade-off between accuracy and compression.
   *
   * **Graph consistency requirement:**  
   * The entire tensor network structure must be specified via @p edges.  
   * The list must contain *all* edges connecting the component tensors of the TPO.  
   * Partial edge specifications (e.g., only nearest neighbors) will lead to undefined behavior.
   *
   * **Ordering requirement:**  
   * The order of tensors in @p O must correspond to the site index order returned by
   * `std::vector<int> Q = GetSiteIndexFromBond(edges)` defined in `tnbp/framework/graph.h`.  
   * That is, for each bond `(site_a, site_b) = edges[m]`, the condition  
   * `Q[i_a] == site_a` ensures that `O[i_a]` represents the TPO on site `site_a`.  
   * This mapping guarantees that each tensor in @p O is correctly associated with its physical site.
   *
   * @tparam TenT Tensor type compatible with TCI operations. Must support contraction, transpose, reshape, and SVD decomposition.
   *
   * @param[in,out] ctx   TCI execution context (including device policy, memory allocator, and backend configuration).
   * @param[in]     edges List of tensor index pairs (i, j) specifying which tensor bonds to process.
   *                      Each pair (i, j) identifies a bond connecting @p O[i] and @p O[j].
   * @param[in,out] O     Tensor product operator, represented as a vector of component tensors.
   *                      This container is updated in place; tensors along processed bonds are replaced with their truncated forms.
   * @param[in]     eps   Singular value cutoff threshold. Singular values ≤ @p eps are discarded.
   *                      Given as a value of type @c real_t<TenT>.
   *
   * @pre
   *  - @p ctx must be properly initialized.
   *  - Each index pair (i, j) in @p edges satisfies 0 ≤ i, j < O.size().
   *  - @p edges must describe the *complete* connectivity of the tensor network.
   *  - Tensor order in @p O must follow the site index sequence obtained from `GetSiteIndexFromBond(edges)`.
   *
   * @post
   *  - @p O is updated in place to reflect reduced bond dimensions.
   *  - A small approximation error may be introduced due to singular value truncation.
   *
   * @note
   *  - The cutoff is “hard”: all singular values ≤ @p eps are discarded.
   *    For relative truncation or maximum-rank-based cutoff, additional logic must be implemented.
   *  - Depending on the backend (CPU/GPU), this function may allocate temporary workspace for SVD.
   *  - The bond-to-site mapping follows `GetSiteIndexFromBond` in `tnbp/framework/graph.h`, and
   *    this must remain consistent with the structure of @p O.
   *
   * @warning
   *  - The operation is in-place; original tensors in @p O will be overwritten.
   *    Make a copy beforehand if the original data must be preserved.
   *  - Incorrect or incomplete edge specifications, or inconsistent ordering of @p O,
   *    may lead to invalid contractions, mismatched ranks, or segmentation faults.
   *
   * @throws std::invalid_argument  If tensor shapes or edge indices are inconsistent.
   * @throws tci::runtime_error     If backend SVD fails or workspace allocation fails.
   *
   * @complexity
   *  For each bond pair, the computational cost is approximately @f$ O(\min(mn^2, m^2 n)) @f$,
   *  where @f$m, n@f$ are the reshaped matrix dimensions. Total cost scales with the number of edge pairs.
   *
   * @see
   *  - tci::svd, tci::contract, tci::reshape
   *  - tnbp::GetSiteIndexFromBond
   *  - MPS/MPO normalization and gate compression routines
   */
  template <typename TenT>
  void OptTPObySVD(context_handle_t<TenT> & ctx,
		   const std::vector<std::pair<int,int>> & edges,
		   std::vector<TenT> & O,
		   real_t<TenT> eps) {
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondIdxT = typename tci::tensor_traits<TenT>::bond_idx_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using SizeT = typename tci::tensor_traits<TenT>::ten_size_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tci::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tci::create_context(ctx_r);

    auto global_sites = GetSiteIndexFromBond(edges);

    for(int edge_address=0; edge_address < edges.size(); edge_address++) {
      int site_a = edges[edge_address].first;
      int site_b = edges[edge_address].second;
      auto it_global_site_address_a = std::find(global_sites.begin(),
						global_sites.end(),
						site_a);
      auto it_global_site_address_b = std::find(global_sites.begin(),
						global_sites.end(),
						site_b);
      auto global_site_address_a = std::distance(global_sites.begin(),
						 it_global_site_address_a);
      auto global_site_address_b = std::distance(global_sites.begin(),
						 it_global_site_address_b);
      std::vector<int> bond_a = GetSurroundingBondIndex(site_a,edges);
      std::vector<int> bond_b = GetSurroundingBondIndex(site_b,edges);
      int bond_address_a;
      int bond_address_b;
      auto it_bond_address_a = std::find(bond_a.begin(),bond_a.end(),
					 edge_address);
      bond_address_a = std::distance(bond_a.begin(),it_bond_address_a);
      auto it_bond_address_b = std::find(bond_b.begin(),bond_b.end(),
					 edge_address);
      bond_address_b = std::distance(bond_b.begin(),it_bond_address_b);

      RankT rank_a = tci::rank(ctx,O[global_site_address_a]);
      RankT rank_b = tci::rank(ctx,O[global_site_address_b]);
      RankT rank_c = rank_a+rank_b-2;
      ShapeT shape_a = tci::shape(ctx,O[global_site_address_a]);
      ShapeT shape_b = tci::shape(ctx,O[global_site_address_b]);
      List<BondLabelT> label_a(rank_a);
      List<BondLabelT> label_b(rank_b);
      List<BondLabelT> label_c(rank_c);
      int uncont_label = 0;
      int k_c = 0;
      for(int k=0; k < rank_a; k++) {
	if( k == bond_address_a ) {
	  label_a[k] = -1;
	} else {
	  label_a[k] = uncont_label;
	  label_c[k_c] = uncont_label;
	  uncont_label++;
	  k_c++;
	}
      }
      for(int k=0; k < rank_b; k++) {
	if( k == bond_address_b ) {
	  label_b[k] = -1;
	} else {
	  label_b[k] = uncont_label;
	  label_c[k_c] = uncont_label;
	  uncont_label++;
	  k_c++;
	}
      }

      TenT C;
      tci::contract(ctx,O[global_site_address_a],label_a,
		    O[global_site_address_b],label_b,
		    C,label_c);
      
      TenT U;
      TenT V;
      RealTenT S;
      TenT D;
      RankT num_rows = rank_a-1;
      RealT trunc_err;
      BondDimT chi_max = shape_a[bond_address_a];
      tci::trunc_svd(ctx,C,num_rows,U,S,V,
		     trunc_err,chi_max,eps);
      tci::for_each(ctx_r,S,[](auto & elem){ elem = std::sqrt(elem); });
      if constexpr (std::is_same_v<TenT,RealTenT>) {
	tci::move(ctx_r,S,D);
      } else {
	D = tci::to_cplx(ctx_r,S);
      }
      tci::diag(ctx,D);
      List<BondLabelT> label_u(rank_a);
      List<BondLabelT> label_d(2);
      for(int k=0; k < rank_a; k++) {
	if( k == bond_address_a ) {
	  label_u[rank_a-1] = -1;
	} else if ( k < bond_address_a ) {
	  label_u[k] = k;
	} else {
	  label_u[k-1] = k;
	}
      }
      label_d[0] = -1;
      label_d[1] = bond_address_a;
      std::iota(label_a.begin(),label_a.end(),0);
      tci::contract(ctx,U,label_u,D,label_d,
		    O[global_site_address_a],label_a);
      List<BondLabelT> label_v(rank_b);
      for(int k=0; k < rank_b; k++) {
	if ( k == bond_address_b ) {
	  label_v[0] = -1;
	} else if( k < bond_address_b ) {
	  label_v[k+1] = k;
	} else {
	  label_v[k] = k;
	}
      }
      label_d[0] = bond_address_b;
      label_d[1] = -1;
      std::iota(label_b.begin(),label_b.end(),0);
      tci::contract(ctx,D,label_d,V,label_v,
		O[global_site_address_b],label_b);
    }
  }
  
}

#endif
