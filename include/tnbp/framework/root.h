/// This file is a part of r-ccs-cms/tnbp
/**
@file tnbp/framework/root.h
@brief Functions to get root of tensor
 */

#ifndef TNBP_FRAMEWORK_ROOT_H
#define TNBP_FRAMEWORK_ROOT_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <utility>

namespace tnbp {

  template <typename TenT>
  void SquareRoot(context_handle_t<TenT> & ctx,
		  const TenT & T,
		  int lb,
		  TenT & S) {

    using BondLabelT  = typename tcapi::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using RealT = typename tcapi::tensor_traits<TenT>::real_t;
    using RealTenT = typename tcapi::tensor_traits<TenT>::real_ten_t;
    using OrderT = typename tcapi::tensor_traits<TenT>::order_t;
    using ShapeT = typename tcapi::tensor_traits<TenT>::shape_t;
    using CoorsT = typename tcapi::tensor_traits<TenT>::elem_coors_t;
    using CtxR = typename tcapi::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tcapi::create_context(ctx_r);

    TenT U;
    TenT V;
    RealTenT E;
    OrderT lb_rt = static_cast<OrderT>(lb);

    tcapi::svd(ctx,T,lb_rt,U,E,V);
    tcapi::for_each(ctx_r,E,[](auto & elem){ elem = std::sqrt(elem); });
    TenT D;
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      D = tcapi::move(ctx_r,E);
    } else {
      D = tcapi::to_cplx(ctx_r,E);
    }

    auto order_U = tcapi::order(ctx,U);
    auto order_V = tcapi::order(ctx,V);
    auto order_D = tcapi::order(ctx,D);
    auto Idx_U = List<BondLabelT>(order_U);
    auto Idx_V = List<BondLabelT>(order_V);
    auto Idx_D = List<BondLabelT>(order_D);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    Idx_U[order_U-1] = -1;
    Idx_D[0] = -1;
    Idx_D[1] = order_U-1;
    auto Idx_O = List<BondLabelT>(order_U);
    std::iota(Idx_O.begin(),Idx_O.end(),0);
    tcapi::contract(ctx,U,Idx_U,D,Idx_D,U,Idx_O);
    std::iota(Idx_U.begin(),Idx_U.end(),0);
    std::iota(Idx_V.begin(),Idx_V.end(),order_U-2);
    Idx_U[order_U-1] = -1;
    Idx_V[0] = -1;
    Idx_O.resize(order_U+order_V-2);
    std::iota(Idx_O.begin(),Idx_O.end(),0);
    tcapi::contract(ctx,U,Idx_U,V,Idx_V,S,Idx_O);

  }

  /**
     General-matrix reference implementation of SquareRootAndInverse based on SVD.

     Kept as the reference path for comparison in tests. Not called by
     SquareRootAndInverse: dense LAPACK SVD (zgesdd/zgesvd) is fragile on the
     Hermitian PSD messenger matrices with wide dynamic range and degenerate
     tails that BP produces (non-termination on some BLAS implementations,
     INFO>0 error returns on others), which is why the production path uses
     the Hermitian eigensolver instead.

     Thresholding (requires sv_min >= 0): singular values sigma <= sv_min are
     dropped from both the root and the inverse root (the inverse-side cut at
     sqrt(sv_min) acts on already-square-rooted values, so the effective
     threshold is the same). For sv_min < 0 the inverse-side comparison is
     against NaN and zeroes S entirely; use the eigh-based
     SquareRootAndInverse, which clamps the floor at zero, instead.
   */
  template <typename TenT>
  void SquareRootAndInverseViaSvd(context_handle_t<TenT> & ctx,
		      const TenT & M,
		      TenT & R,
		      TenT & S,
		      real_t<TenT> sv_min) {

    using BondLabelT  = typename tcapi::tensor_traits<TenT>::bond_label_t;
    using RealTenT = typename tcapi::tensor_traits<TenT>::real_ten_t;
    using OrderT = typename tcapi::tensor_traits<TenT>::order_t;
    using CtxR = typename tcapi::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tcapi::create_context(ctx_r);

    TenT U;
    TenT V;
    RealTenT E;
    OrderT lb_rt = static_cast<OrderT>(1);
    tcapi::svd(ctx,M,lb_rt,U,E,V);
    TenT D;
    TenT F;
    tcapi::for_each(ctx_r,E,[&sv_min](auto & elem){
      if( std::abs(elem) > sv_min ) { elem = std::sqrt(elem); }
      else { elem = 0.0; } });
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      D = tcapi::copy(ctx_r,E);
    } else {
      D = tcapi::to_cplx(ctx_r,E);
    }
    tcapi::for_each(ctx_r,E,[&sv_min](auto & elem) {
      if( std::abs(elem) > std::sqrt(sv_min)) { elem = elem/(elem*elem); }
      else { elem = 0.0; } });
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      F = tcapi::copy(ctx_r,E);
    } else {
      F = tcapi::to_cplx(ctx_r,E);
    }

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
    tcapi::contract(ctx,U,Idx_U,D,Idx_D,W,Idx_R);
    Idx_U[0] = 0;
    Idx_U[1] = -1;
    Idx_V[0] = -1;
    Idx_V[1] = 1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    tcapi::contract(ctx,W,Idx_U,V,Idx_D,R,Idx_R);

    tcapi::cplx_conj(ctx,U);
    tcapi::cplx_conj(ctx,V);
    Idx_V[0] = -1;
    Idx_V[1] = 0;
    Idx_D[0] = 1;
    Idx_D[1] = -1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    tcapi::contract(ctx,F,Idx_D,V,Idx_V,W,Idx_R);
    Idx_V[0] = 0;
    Idx_V[1] = -1;
    Idx_U[0] = 1;
    Idx_U[1] = -1;
    Idx_R[0] = 0;
    Idx_R[1] = 1;
    tcapi::contract(ctx,W,Idx_V,U,Idx_U,S,Idx_R);
  }

  /**
     Build a diagonal 2-index tensor by mapping func over a copy of the real
     diagonal spectrum tensor E, shared by the root / inverse-root constructions so the
     threshold mapping exists once.
   */
  template <typename TenT, typename Func>
  TenT DiagFromSpectrum(context_handle_t<TenT> & ctx,
			context_handle_t<real_ten_t<TenT>> & ctx_r,
			const real_ten_t<TenT> & E,
			Func && func) {
    using RealTenT = typename tcapi::tensor_traits<TenT>::real_ten_t;
    RealTenT Emapped = tcapi::copy(ctx_r,E);
    tcapi::for_each_with_coors(ctx_r,Emapped,[&](auto& elem, const auto& coors) {
      if (coors[0] == coors[1]) func(elem);
    });
    TenT D;
    if constexpr (std::is_same_v<TenT,RealTenT>) {
      // Emapped is a private local, so ownership can be transferred directly.
      D = std::move(Emapped);
    } else {
      D = tcapi::to_cplx(ctx_r,Emapped);
    }
    return D;
  }

  /**
     Compute the matrix square root R = M^{1/2} and the thresholded inverse
     square root S = M^{-1/2} of a 2-index tensor M.

     Precondition: M must be Hermitian positive-semidefinite up to
     floating-point rounding. Every in-library producer of BP messenger
     tensors satisfies this by construction (identity initialization,
     conjugate-paired sandwich updates, nonnegative-diagonal resets), so a
     violation indicates a caller bug. The precondition is not checked at
     runtime, matching the guard level of the SVD-based reference
     implementation: rounding-level non-Hermitian components are absorbed by
     the symmetrization below, and negative eigenvalues by the clamping.

     The decomposition uses the Hermitian eigensolver (tcapi::eigh) instead of
     general SVD: dense LAPACK SVD is unreliable on Hermitian PSD matrices
     whose spectrum spans many decades with a degenerate tail, while the
     symmetric eigensolver family handles them robustly. Rounding-level
     non-Hermitian components are removed by symmetrizing M before the
     decomposition.

     Eigenvalues lambda <= max(sv_min, 0) are dropped from both R and S.
     The signed comparison also clamps numerically negative eigenvalues of a
     PSD-up-to-rounding input to zero, keeping sqrt away from negative
     arguments. For exactly PSD input and sv_min >= 0 this matches the
     SVD-based reference implementation (SquareRootAndInverseViaSvd), since
     R and S are matrix functions of M and therefore independent of the
     factorization choice.
   */
  template <typename TenT>
  void SquareRootAndInverse(context_handle_t<TenT> & ctx,
		      const TenT & M,
		      TenT & R,
		      TenT & S,
		      real_t<TenT> sv_min) {

    using BondLabelT  = typename tcapi::tensor_traits<TenT>::bond_label_t;
    using ElemT = typename tcapi::tensor_traits<TenT>::elem_t;
    using RealT = typename tcapi::tensor_traits<TenT>::real_t;
    using RealTenT = typename tcapi::tensor_traits<TenT>::real_ten_t;
    using OrderT = typename tcapi::tensor_traits<TenT>::order_t;
    using CtxR = typename tcapi::tensor_traits<RealTenT>::context_handle_t;
    CtxR ctx_r;
    tcapi::create_context(ctx_r);

    TenT Mt = tcapi::copy(ctx,M);
    tcapi::transpose(ctx,Mt,{1,0});
    tcapi::cplx_conj(ctx,Mt);

    // Symmetrize to remove rounding-level non-Hermitian components before
    // handing the matrix to the Hermitian eigensolver.
    TenT Msym = tcapi::linear_combine<TenT>(ctx,{std::cref(M),std::cref(Mt)},
					  {ElemT(0.5),ElemT(0.5)});

    RealTenT E;
    TenT V;
    OrderT lb_rt = static_cast<OrderT>(1);
    tcapi::eigh(ctx,Msym,lb_rt,E,V);

    // Signed threshold: eigenvalues at or below the clamped floor are
    // dropped, which also projects numerically negative eigenvalues to zero.
    // The drop rule exists once so R and S can never disagree on the cut.
    RealT floor_ev = std::max(sv_min,RealT(0.0));
    auto thresholded_diag = [&](auto && map) {
      return DiagFromSpectrum<TenT>(ctx,ctx_r,E,[&](auto & elem){
	if( elem > floor_ev ) { elem = map(elem); }
	else { elem = 0.0; } });
    };
    TenT D = thresholded_diag([](auto x){ return std::sqrt(x); });
    TenT F = thresholded_diag([](auto x){ return 1.0/std::sqrt(x); });

    TenT Vc = tcapi::copy(ctx,V);
    tcapi::cplx_conj(ctx,Vc);

    // R = V diag(sqrt(lambda)) V^dagger, S = V diag(1/sqrt(lambda)) V^dagger.
    // The last bond of V enumerates eigenvectors, so the adjoint orientation
    // of the conjugated factor is supplied by the contraction labels alone.
    auto sandwich_with_v = [&](const TenT & Dg, TenT & Out) {
      List<BondLabelT> Idx_L = {0, -1};
      List<BondLabelT> Idx_Dg = {-1, 1};
      List<BondLabelT> Idx_Vc = {1, -1};
      List<BondLabelT> Idx_O = {0, 1};
      TenT W;
      tcapi::contract(ctx,V,Idx_L,Dg,Idx_Dg,W,Idx_O);
      tcapi::contract(ctx,W,Idx_L,Vc,Idx_Vc,Out,Idx_O);
    };
    sandwich_with_v(D,R);
    sandwich_with_v(F,S);


  }


}

#endif
