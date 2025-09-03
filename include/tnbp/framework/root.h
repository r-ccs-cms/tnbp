/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/root.h
@brief Functions to get root of tensor
 */

#ifndef TNBP_FRAMEWORK_ROOT_H
#define TNBP_FRAMEOWRK_ROOT_H

namespace tnbp {

  template <typename TenT>
  void SquareRoot(tci::context_handle_t<TenT> & ctx,
		  const TenT & T,
		  int lb,
		  TenT & S) {

    using IntT  = typename tci::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealT = typename tci::tensor_traits<TenT>::real_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;

    TenT U;
    TenT V;
    RealTenT D;
    RankT lb_rt = static_cast<RankT>(lb);

    tci::svd(ctx,T,lb_rt,U,D,V);
    tci::for_each(ctx,D,[](ElemT & elem) { elem = std::sqrt(elem); });
    tci::diag(ctx,D);

    auto shape_U = tci::shape(ctx,U);
    auto shape_V = tci::shape(ctx,V);
    auto shape_D = tci::shape(ctx,D);
    auto Idx_U = shape_U;
    auto Idx_V = shape_V;
    auto Idx_D = shape_D;
    auto Rank_U = tci::rank(ctx,U);
    auto Rank_V = tci::rank(ctx,V);
    auto Rank_D = tci::rank(ctx,D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    tci::contract(ctx,U,D,Idx_U,Idx_D,U);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),Rank_U-1);
    Idx_U[Rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,U,V,Idx_U,Idx_D,S);
    
  }

  template <typename TenT, typename RealT>
  void SquareRootAndInverse(const TenT & M,
		      TenT & R,
		      TenT & S,
		      RealT sv_min) {
    
    using IntT  = typename tci::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using RealTenT = typename tci::tensor_traits<TenT>::real_ten_t;
    using RankT = typename tci::tensor_traits<TenT>::rank_t;

    TenT U;
    TenT V;
    RealTenT D;
    RankT lb_rt = static_cast<RankT>(lb);

    tci::svd(ctx,T,lb_rt,U,D,V);
    tci::for_each(ctx,D,[](ElemT & elem) { elem = std::sqrt(elem); });
    auto F = D;
    tci::diag(ctx,D);
    tci::for_each(ctx,F,[](ElemT & elem) {
      if( std::abs(elem) > sv_min) elem = 1.0/elem;
      else elem = ElemT(0.0); });
    tci::diag(ctx,F);
    
    auto shape_U = tci::shape(ctx,U);
    auto shape_V = tci::shape(ctx,V);
    auto shape_D = tci::shape(ctx,D);
    auto Idx_U = shape_U;
    auto Idx_V = shape_V;
    auto Idx_D = shape_D;
    auto Rank_U = tci::rank(ctx,U);
    auto Rank_V = tci::rank(ctx,V);
    auto Rank_D = tci::rank(ctx,D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    auto W = U;
    tci::contract(ctx,W,D,Idx_U,Idx_D,W);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),Rank_U-1);
    Idx_U[Rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,W,V,Idx_U,Idx_D,R);

    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    tci::contract(ctx,U,F,Idx_U,Idx_D,U);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),Rank_U-1);
    Idx_U[Rank_U-1] = -1;
    Idx_V[0] = -1;
    tci::contract(ctx,U,V,Idx_U,Idx_D,S);
  }
		    
  
}

#endif
