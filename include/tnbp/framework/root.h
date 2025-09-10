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

    TenT U;
    TenT V;
    RealTenT S;
    RankT lb_rt = static_cast<RankT>(lb);

    tci::svd(ctx,T,lb_rt,U,S,V);
    ShapeT shape_D = tci::shape(ctx,S);
    std::vector<ElemT> data_D(shape_D[0]);
    auto it_D = data_D.begin();
    tci::for_each(ctx,S,[&it_D](RealT & elem) {
      *it_D++ = static_cast<ElemT>(std::sqrt(elem)); });
    TenT D;
    tci::allocate(ctx,shape_D,D);
    it_D = data_D.begin();
    tci::for_each(ctx,D,[&it_D](ElemT & elem) {
      elem = *it_D++; });

    auto rank_U = tci::rank(ctx,U);
    auto rank_V = tci::rank(ctx,V);
    auto rank_D = tci::rank(ctx,D);
    auto Idx_U = std::vector<BondLabelT>(rank_U);
    auto Idx_V = std::vector<BondLabelT>(rank_V);
    auto Idx_D = std::vector<IntT>(rank_D);
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

    TenT U;
    TenT V;
    RealTenT S;
    RankT lb_rt = static_cast<RankT>(lb);
    tci::svd(ctx,T,lb_rt,U,D,V);
    ShapeT shape_S = tci::shape(ctx,S);
    std::vector<ElemT> data_D(shape_S[0]);
    std::vector<ElemT> data_F(shape_S[0]);
    auto it_D = data_D.begin();
    tci::for_each(ctx,S,[&it_D,sv_min](const RealT & elem) {
      if( std::abs(elem) > sv_min ) {
	*it_D++ = std::sqrt(elem);
      } else {
	*it_D++ = static:cast<ElemT>(0.0);
      } });
    auto it_F = data_F.begin();
    tci::for_each(ctx,S,[&it_F,sv_min](const RealT & elem) {
      if( std::abs(elem) > sv_min ) {
	*it_F++ = static_cast<ElemT>(1.0/elem);
      } else {
	*it_F++ = static_cast<ElemT>(0.0);
      } });
    TenT D; tci::allocate(ctx,shape_S,D);
    TenT F; tci::allocate(ctx,shape_S,F);
    tci::diag(ctx,D);
    tci::diag(ctx,F);
    
    auto Rank_U = tci::rank(ctx,U);
    auto Rank_V = tci::rank(ctx,V);
    auto Rank_D = tci::rank(ctx,D);
    auto Idx_U = std::vector<BondLabelT>(rank_U);
    auto Idx_V = std::vector<BondLabelT>(rank_V);
    auto Idx_D = std::vector<BondLabelT>(rank_D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[Rank_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = Rank_U-1;
    auto W = U;
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
