/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/root.h
@brief Functions to get root of tensor
 */

#ifndef TNBP_FRAMEWORK_ROOT_H
#define TNBP_FRAMEWORK_ROOT_H

#include <type_traits>

namespace tnbp {

  template <typename TenT>
  void SquareRoot(context_handle_t<TenT> & ctx,
		  const TenT & T,
		  int lb,
		  TenT & S) {

    using BondLabelT  = typename tci::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tci::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tci::create_context(ctx_r);

    TenT U;
    TenT V;
    RealTenT E;
    RankT lb_rt = static_cast<RankT>(lb);

    tci::svd(ctx,T,lb_rt,U,E,V);
    tci::for_each(ctx_r,E,[](auto & elem){ elem = std::sqrt(elem); });
    TenT D;
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      tci::move(ctx_r,E,D);
    } else {
      tci::to_cplx(ctx_r,E,D);
    }
    tci::diag(ctx,D);

    auto rank_U = tci::rank(ctx,U);
    auto rank_V = tci::rank(ctx,V);
    auto rank_D = tci::rank(ctx,D);
    auto Idx_U = List<BondLabelT>(rank_U);
    auto Idx_V = List<BondLabelT>(rank_V);    
    auto Idx_D = List<BondLabelT>(rank_D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = rank_U-1;
    tci::contract(ctx,U,Idx_U,D,Idx_D,U);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),rank_U-1);
    Idx_U[rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,U,Idx_U,V,Idx_D,S);
    
  }

  template <typename TenT>
  void SquareRootAndInverse(context_handle_t<TenT> & ctx,
		      const TenT & M,
		      TenT & R,
		      TenT & S,
		      real_t<TenT> sv_min) {
    
    using BondLabelT  = typename tci::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tci::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tci::create_context(ctx_r);

    TenT U;
    TenT V;
    RealTenT E;
    RankT lb_rt = static_cast<RankT>(1);
    tci::svd(ctx,M,lb_rt,U,E,V);
    TenT D;
    TenT F;
    tci::for_each(ctx_r,E,[&sv_min](auto & elem){
      if( std::abs(elem) > sv_min ) { elem = std::sqrt(elem); }
      else { elem = 0.0; } });
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      tci::copy(ctx_r,E,D);
    } else {
      tci::to_cplx(ctx_r,E,D);
    }
    tci::for_each(ctx_r,E,[&sv_min](auto & elem) {
      if( std::abs(elem) > std::sqrt(sv_min)) { elem = elem/(elem*elem); }
      else { elem = 0.0; } });
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      tci::copy(ctx_r,E,F);
    } else {
      tci::to_cplx(ctx_r,E,F);
    }
    tci::diag(ctx,D);
    tci::diag(ctx,F);
    
    auto Rank_U = tci::rank(ctx,U);
    auto Rank_V = tci::rank(ctx,V);
    auto Rank_D = tci::rank(ctx,D);
    auto Idx_U = List<BondLabelT>(2);
    auto Idx_V = List<BondLabelT>(2);
    auto Idx_D = List<BondLabelT>(2);
    auto Idx_R = List<BondLabelT>(2);
    Idx_U[0] = 0;
    Idx_U[1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = 1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    TenT W;
    tci::contract(ctx,U,Idx_U,D,Idx_D,W,Idx_R);
    Idx_U[0] = 0;
    Idx_U[1] = -1;
    Idx_V[0] = -1;
    Idx_V[1] = 1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    tci::contract(ctx,W,Idx_U,V,Idx_D,R,Idx_R);

    tci::cplx_conj(ctx,U);
    tci::cplx_conj(ctx,V);
    Idx_V[0] = -1;
    Idx_V[1] = 0;
    Idx_D[0] = 1;
    Idx_D[1] = -1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    tci::contract(ctx,F,Idx_D,V,Idx_V,W,Idx_R);
    Idx_V[0] = 0;
    Idx_V[1] = -1;
    Idx_U[0] = 1;
    Idx_U[1] = -1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    tci::contract(ctx,W,Idx_V,U,Idx_U,S,Idx_R);
  }
		    
  
}

#endif
