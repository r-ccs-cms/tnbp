/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/root.h
@brief Functions to get root of tensor
 */

#ifndef TNBP_FRAMEWORK_ROOT_H
#define TNBP_FRAMEOWRK_ROOT_H

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
    tci::create_context<RealTenT>(ctx_r);

    TenT U;
    TenT V;
    RealTenT S;
    RankT lb_rt = static_cast<RankT>(lb);

    tci::svd(ctx,T,lb_rt,U,S,V);
    tci::for_each(ctx_r,S,[](RealT & elem){ elem = std::sqrt(elem); });
    TenT D;
    tci::convert(ctx_r,S,ctx,D);
    tci::diag(ctx,D);

    auto rank_U = tci::rank(ctx,U);
    auto rank_V = tci::rank(ctx,V);
    auto rank_D = tci::rank(ctx,D);
    auto Idx_U = List<BondLabelT>(rank_U);
    auto Idx_V = List<BondLabelT>(rank_V);    
    auto Idx_D = List<BondLabelT>(rank_D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    tci::contract(ctx,U,Idx_U,D,Idx_D,U);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),Rank_U-1);
    Idx_U[Rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,U,Idx_U,V,Idx_D,S);
    
  }

  template <typename TenT>
  void SquareRootAndInverse(const TenT & M,
		      TenT & R,
		      TenT & S,
		      real_t<TenT> sv_min) {
    
    using BondLabelT  = typename tci::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tci::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tci::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctr_r;
    tci::create_context<RealTenT>(ctx_r);

    TenT U;
    TenT V;
    RealTenT S;
    RankT lb_rt = static_cast<RankT>(lb);
    tci::svd(ctx,T,lb_rt,U,D,V);
    tci::for_each(ctx_r,S,[](RealT & elem){
      if( elem > sv_min ) { elem = std::sqrt(elem); }
      else { elem = 0.0; } });
    TenT D;
    tci::convert(ctx_r,S,ctx,D);
    tci::for_each(ctx_r,S,[](RealT & elem) {
      if( elem > std::sqrt(sv_min)) { elem = 1.0/elem; }
      else { elem = 0.0; } });
    TenT F;
    tci::convert(ctx_r,S,ctx,F);
    tci::diag(ctx,D);
    tci::diag(ctx,F);
    
    auto Rank_U = tci::rank(ctx,U);
    auto Rank_V = tci::rank(ctx,V);
    auto Rank_D = tci::rank(ctx,D);
    auto Idx_U = List<BondLabelT>(rank_U);
    auto Idx_V = List<BondLabelT>(rank_V);
    auto Idx_D = List<BondLabelT>(rank_D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    TenT W;
    tci::copy(ctx,U,W);
    tci::contract(ctx,W,Idx_U,D,Idx_D,W);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),Rank_U-1);
    Idx_U[Rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,W,Idx_U,V,Idx_D,R);

    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    tci::contract(ctx,U,Idx_U,F,Idx_D,U);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),Rank_U-1);
    Idx_U[Rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,U,Idx_U,V,Idx_D,S);
  }
		    
  
}

#endif
