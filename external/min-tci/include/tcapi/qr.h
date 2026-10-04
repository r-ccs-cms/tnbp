#pragma once
#include "tensor_manipulation.h"
#include "miscellaneous.h"
#include "gqten/tensor/out_of_place_ops/qr.h"

namespace tcapi {
namespace detail {
template<class TenT>
void validate_qr(const TenT& a, order_t<TenT> rows, const TenT& left, const TenT& right) {
    if (rows<1 || static_cast<std::size_t>(rows)>=a.Rank())
        throw std::invalid_argument("min-tcapi: QR/LQ requires 1 <= row bonds < rank");
    if (std::addressof(left)==std::addressof(right))
        throw std::invalid_argument("min-tcapi: QR/LQ outputs must be distinct");
    const auto count=checked_size<TenT>(a.Shape());
    if (count>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("min-tcapi: QR/LQ exceeds gqten integer size limit");
    for (std::size_t i=0;i<count;++i) {
        const auto x=numeric_value(a,i);
        if (!std::isfinite(std::real(x)) || !std::isfinite(std::imag(x)))
            throw std::domain_error("min-tcapi: QR/LQ requires finite input");
    }
}
} // namespace detail

template<class TenT>
void qr(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows, TenT& q, TenT& r) {
    detail::verbose::call diagnostic("qr", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    detail::validate_qr(a,rows,q,r);
    TenT new_q,new_r;
    gqten::QR(&a,static_cast<std::size_t>(rows),&new_q,&new_r);
    detail::replace(q,std::move(new_q));
    detail::replace(r,std::move(new_r));
}

template<class TenT>
void lq(context_handle_t<TenT>& ctx, const TenT& a, order_t<TenT> rows, TenT& l, TenT& q) {
    detail::verbose::call diagnostic("lq", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    detail::validate_qr(a,rows,l,q);
    // A^T = Q' R' implies A = R'^T Q'^T, including complex tensors.
    // Move column axes before row axes while preserving order within each group.
    const auto cols=a.Rank()-static_cast<std::size_t>(rows);
    List<bond_idx_t<TenT>> permutation;
    for (std::size_t i=rows;i<a.Rank();++i) permutation.push_back(i);
    for (order_t<TenT> i=0;i<rows;++i) permutation.push_back(i);
    auto transposed=tcapi::copy(ctx,a);
    transposed.Transpose(permutation);
    TenT new_q,new_l;
    gqten::QR(&transposed,cols,&new_q,&new_l);
    permutation.clear();
    for (order_t<TenT> i=1;i<=rows;++i) permutation.push_back(i);
    permutation.push_back(0);
    new_l.Transpose(permutation);
    permutation.clear();
    permutation.push_back(cols);
    for (std::size_t i=0;i<cols;++i) permutation.push_back(i);
    new_q.Transpose(permutation);
    detail::replace(l,std::move(new_l));
    detail::replace(q,std::move(new_q));
}
} // namespace tcapi
