#pragma once
#include "miscellaneous.h"
#include "gqten/tensor/out_of_place_ops/svd.h"
#include <algorithm>

namespace tcapi {
namespace detail {
template<class TenT>
std::size_t validate_svd(const TenT& a, order_t<TenT> rows,
                         TenT& u, real_ten_t<TenT>& sigma, TenT& vd) {
    if (rows<1 || static_cast<std::size_t>(rows)>=a.Rank())
        throw std::invalid_argument("min-tcapi: SVD requires 1 <= row bonds < rank");
    const void* up=std::addressof(u); const void* sp=std::addressof(sigma);
    const void* vp=std::addressof(vd);
    if (up==sp || up==vp || sp==vp)
        throw std::invalid_argument("min-tcapi: SVD outputs must be distinct");
    const auto count=checked_size<TenT>(a.Shape());
    if (count>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("min-tcapi: SVD exceeds gqten integer size limit");
    std::size_t m=1;
    for (order_t<TenT> i=0;i<rows;++i) m*=a.Shape()[i];
    const auto k=std::min(m,count/m);
    checked_size<real_ten_t<TenT>>({static_cast<bond_dim_t<TenT>>(k),static_cast<bond_dim_t<TenT>>(k)});
    for (std::size_t i=0;i<count;++i) {
        const auto x=numeric_value(a,i);
        if (!std::isfinite(std::real(x)) || !std::isfinite(std::imag(x)))
            throw std::domain_error("min-tcapi: SVD requires finite input");
    }
    return k;
}
template<class TenT>
auto full_svd(const TenT& a, order_t<TenT> rows, std::size_t k, TenT& u, TenT& vd) {
    using R=real_t<TenT>;
    R* raw=nullptr; std::size_t actual=0; R unused_error{};
    // Unlike SVD, TruncSVD reports LAPACK failure. Keeping k values disables
    // truncation. Explicit gesvd avoids the native auto-driver retry copy bug.
    const auto info=gqten::TruncSVD(&a,static_cast<std::size_t>(rows),R(0),k,k,
        &u,&vd,raw,&actual,&unused_error,R(0),"svd");
    std::unique_ptr<R,decltype(&std::free)> values(raw,&std::free);
    if (info!=0) throw std::runtime_error("min-tcapi: gqten SVD failed");
    if (actual!=k) throw std::runtime_error("min-tcapi: unexpected SVD dimension");
    for (std::size_t i=0;i<k;++i)
        if (!std::isfinite(values.get()[i]) || values.get()[i]<R(0))
            throw std::runtime_error("min-tcapi: invalid singular values");
    return values;
}
template<class TenT>
real_ten_t<TenT> singular_diagonal(const real_t<TenT>* values, std::size_t k) {
    const auto dim=static_cast<bond_dim_t<TenT>>(k);
    return construct<real_ten_t<TenT>>({dim,dim},[&](std::size_t i) {
        return i%k==i/k ? values[i%k] : real_t<TenT>(0);
    });
}
} // namespace detail

template<class TenT>
void svd(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows,
         TenT& u, real_ten_t<TenT>& sigma, TenT& v_dag) {
    detail::verbose::call diagnostic("svd", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    const auto k=detail::validate_svd(a,rows,u,sigma,v_dag);
    TenT new_u,new_v;
    auto values=detail::full_svd(a,rows,k,new_u,new_v);
    auto new_s=detail::singular_diagonal<TenT>(values.get(),k);
    detail::replace(u,std::move(new_u));
    detail::replace(sigma,std::move(new_s));
    detail::replace(v_dag,std::move(new_v));
}

template<class TenT>
void trunc_svd(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows,
               TenT& u, real_ten_t<TenT>& sigma, TenT& v_dag, real_t<TenT>& trunc_err,
               bond_dim_t<TenT> chi_min, bond_dim_t<TenT> chi_max,
               real_t<TenT> target, real_t<TenT> s_min) {
    detail::verbose::call diagnostic("trunc_svd", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
        detail::verbose::field(log, "chi_min", chi_min);
        detail::verbose::field(log, "chi_max", chi_max);
        detail::verbose::field(log, "target", target);
        detail::verbose::field(log, "s_min", s_min);
    });
    using R=real_t<TenT>;
    detail::require_context(ctx);
    if (chi_min<1 || chi_max<chi_min || !std::isfinite(target) || target<R(0)
        || !std::isfinite(s_min) || s_min<R(0))
        throw std::invalid_argument("min-tcapi: invalid SVD truncation parameters");
    const auto k=detail::validate_svd(a,rows,u,sigma,v_dag);
    TenT full_u,full_v;
    auto values=detail::full_svd(a,rows,k,full_u,full_v);
    std::size_t survivors=k;
    while (survivors && values.get()[survivors-1]<s_min) --survivors;
    if (!survivors)
        throw std::domain_error("min-tcapi: truncation to dimension zero is unsupported");
    // Normalize before squaring; suffix sums avoid cancellation when the
    // discarded weight is much smaller than the retained weight.
    std::vector<long double> tail(k+1,0);
    const long double largest=values.get()[0];
    for (std::size_t i=k;i>0;--i) {
        const long double ratio=largest==0 ? 0 : values.get()[i-1]/largest;
        tail[i-1]=tail[i]+ratio*ratio;
    }
    auto error=[&](std::size_t chi) { return tail[0]==0 ? 0.L : tail[chi]/tail[0]; };
    std::size_t chi=std::min(survivors,static_cast<std::size_t>(chi_min));
    const auto cap=std::min(survivors,static_cast<std::size_t>(chi_max));
    while (chi<cap && error(chi)>static_cast<long double>(target)) ++chi;
    const R new_error=static_cast<R>(error(chi));
    auto us=full_u.Shape(),vs=full_v.Shape();
    us.back()=static_cast<bond_dim_t<TenT>>(chi);
    vs.front()=static_cast<bond_dim_t<TenT>>(chi);
    auto new_u=detail::construct<TenT>(us,[&](std::size_t i) { return detail::numeric_value(full_u,i); });
    auto new_v=detail::construct<TenT>(vs,[&](std::size_t i) {
        return detail::numeric_value(full_v,i%chi+(i/chi)*k);
    });
    auto new_s=detail::singular_diagonal<TenT>(values.get(),chi);
    detail::replace(u,std::move(new_u));
    detail::replace(sigma,std::move(new_s));
    detail::replace(v_dag,std::move(new_v));
    trunc_err=new_error;
}

template<class TenT>
void trunc_svd(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows,
               TenT& u, real_ten_t<TenT>& sigma, TenT& v_dag, real_t<TenT>& trunc_err,
               bond_dim_t<TenT> chi_max, real_t<TenT> s_min) {
    detail::verbose::call diagnostic("trunc_svd", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
        detail::verbose::field(log, "chi_max", chi_max);
        detail::verbose::field(log, "s_min", s_min);
    });
    tcapi::trunc_svd(ctx,a,rows,u,sigma,v_dag,trunc_err,1,chi_max,real_t<TenT>(0),s_min);
}
} // namespace tcapi
