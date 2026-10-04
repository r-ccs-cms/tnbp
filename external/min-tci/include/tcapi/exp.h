#pragma once
#include "eigh.h"
#include "gqten/tensor/out_of_place_ops/exp.h"

namespace tcapi {
// This backend implements the symmetric/Hermitian matrix exponential only.
template<class TenT>
void exp(context_handle_t<TenT>& ctx, const TenT& in, order_t<TenT> rows, TenT& out) {
    detail::verbose::call diagnostic("exp", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    detail::validate_eigh(in,rows);
    TenT result;
    // Native implementation forms B=V exp(lambda/2), then B B^dagger.
    gqten::ExpHermExact(&in,static_cast<std::size_t>(rows),&result);
    detail::replace(out,std::move(result));
}
template<class TenT>
void exp(context_handle_t<TenT>& ctx, TenT& inout, order_t<TenT> rows) {
    detail::verbose::call diagnostic("exp", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "rows", rows);
    });
    tcapi::exp(ctx,static_cast<const TenT&>(inout),rows,inout);
}
} // namespace tcapi
