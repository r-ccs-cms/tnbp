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
  void Apply(context_handle_t<TenT> & ctx,
	     std::vector<TenT> & V,
	     const std::vector<int> & SiteIdx,
	     const std::vector<TenT> & O) {
    
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using BondLabelT = typename tci::tensor_traits<TenT>::bond_label_t;
    using BondDimT = typename tci::tensor_traits<TenT>::bond_dim_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;

    for(int address=0; address < SiteIdx.size(); address++) {
      int site = SiteIdx[address];
      TenT W;
      RankT rank_v = tci::rank(ctx,V[address]);
      RankT rank_o = tci::rank(ctx,O[site]);
      RankT num_vb = rank_v-1;
      std::vector<BondLabelT> label_v(rank_v);
      std::vector<BondLabelT> label_o(rank_o);
      std::vector<BondLabelT> label_r(rank_v+rank_o-2);
      ShapeT new_shape_v(rank_v);
      ShapeT shape_v = tci::shape(ctx,V[address]);
      ShapeT shape_o = tci::shape(ctx,O[site]);
      for(RankT k=0; k < num_vb; k++) {
	label_v[k] = 2*k+0;
	label_o[k] = 2*k+1;
	label_r[2*k+0] = 2*k+0;
	label_r[2*k+1] = 2*k+1;
	new_shape_v[k] = shape_v[k]*shape_o[k];
      }
      label_v[num_vb] = -1;
      label_o[num_vb] = 2*num_vb+1;
      label_o[num_vb+1] = -1;
      label_r[2*num_vb] = 2*num_vb+1;
      new_shape_v[k] = shape_o[num_vb];
      tci::contract(ctx,O[site],label_o,V[address],label_v,V[address],label_r);
      tci::reshape(ctx,V[address],new_shape_v);
    }
  }
}

#endif

