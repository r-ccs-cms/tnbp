#pragma once
#include "miscellaneous.h"
#include "gqten/tensor/out_of_place_ops/eig.h"

namespace tcapi {
namespace detail {
template<class TenT>
std::size_t validate_eigh(const TenT& a, order_t<TenT> rows) {
    if (rows<1 || static_cast<std::size_t>(rows)>=a.Rank())
        throw std::invalid_argument("min-tcapi: Hermitian eigensolver requires 1 <= row bonds < rank");
    const auto count=checked_size<TenT>(a.Shape());
    if (count>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("min-tcapi: Hermitian eigensolver exceeds gqten integer size limit");
    std::size_t n=1;
    for (order_t<TenT> i=0;i<rows;++i) n*=a.Shape()[i];
    if (n!=count/n)
        throw std::invalid_argument("min-tcapi: Hermitian eigensolver requires a square matricization");
    for (std::size_t i=0;i<count;++i) {
        const auto x=numeric_value(a,i);
        if (!std::isfinite(std::real(x)) || !std::isfinite(std::imag(x)))
            throw std::domain_error("min-tcapi: Hermitian eigensolver requires finite input");
    }
    return n;
}
template<class TenT>
auto hermitian_factors(const TenT& a, order_t<TenT> rows, std::size_t n, bool vectors) {
    // Materialize logical values: native EigHerm applies only real(scale)
    // after diagonalization, which also reverses ordering for negative scale.
    auto materialized=construct<TenT>(a.Shape(),[&](std::size_t i){return numeric_value(a,i);});
    real_t<TenT>* w=nullptr; elem_t<TenT>* v=nullptr; std::size_t actual=0;
    gqten::EigHerm(&materialized,static_cast<std::size_t>(rows),w,v,&actual,vectors?'V':'N','U');
    std::unique_ptr<real_t<TenT>,decltype(&std::free)> values(w,&std::free);
    std::unique_ptr<elem_t<TenT>,decltype(&std::free)> eigenvectors(v,&std::free);
    if (actual!=n) throw std::runtime_error("min-tcapi: unexpected eigenvalue count");
    return std::make_pair(std::move(values),std::move(eigenvectors));
}
} // namespace detail

template<class TenT>
void eigh(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows,
          real_ten_t<TenT>& lambda_mat, TenT& v) {
    detail::verbose::call diagnostic("eigh", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    if (static_cast<const void*>(std::addressof(lambda_mat))==static_cast<const void*>(std::addressof(v)))
        throw std::invalid_argument("min-tcapi: eigh outputs must be distinct");
    const auto n=detail::validate_eigh(a,rows);
    auto factors=detail::hermitian_factors(a,rows,n,true);
    const auto dim=static_cast<bond_dim_t<TenT>>(n);
    auto new_w=detail::construct<real_ten_t<TenT>>({dim,dim},[&](std::size_t i) {
        return i%n==i/n ? factors.first.get()[i%n] : real_t<TenT>(0);
    });
    shape_t<TenT> vs(a.Shape().begin(),a.Shape().begin()+rows);vs.push_back(dim);
    auto new_v=detail::construct<TenT>(vs,[&](std::size_t i){return factors.second.get()[i];});
    detail::replace(lambda_mat,std::move(new_w));
    detail::replace(v,std::move(new_v));
}

template<class TenT>
void eigvalsh(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows,
              real_ten_t<TenT>& w) {
    detail::verbose::call diagnostic("eigvalsh", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    const auto n=detail::validate_eigh(a,rows);
    auto factors=detail::hermitian_factors(a,rows,n,false);
    auto new_w=detail::construct<real_ten_t<TenT>>({static_cast<bond_dim_t<TenT>>(n)},
        [&](std::size_t i){return factors.first.get()[i];});
    detail::replace(w,std::move(new_w));
}
} // namespace tcapi
