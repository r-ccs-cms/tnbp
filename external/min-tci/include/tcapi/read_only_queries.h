#pragma once
#include "detail/core.h"
namespace tcapi {
template<class TenT> order_t<TenT> order(context_handle_t<TenT>& ctx, const TenT& a) {
    detail::verbose::call diagnostic("order", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(a.Shape());
    return static_cast<order_t<TenT>>(a.Rank());
}
template<class TenT> shape_t<TenT> shape(context_handle_t<TenT>& ctx, const TenT& a) {
    detail::verbose::call diagnostic("shape", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    detail::require_context(ctx);
    return a.Shape();
}
template<class TenT> ten_size_t<TenT> size(context_handle_t<TenT>& ctx, const TenT& a) {
    detail::verbose::call diagnostic("size", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    detail::require_context(ctx);
    return detail::checked_size<TenT>(a.Shape());
}
template<class TenT> std::size_t size_bytes(context_handle_t<TenT>& ctx, const TenT& a) {
    detail::verbose::call diagnostic("size_bytes", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    // Numeric payload only: the scalar lives inline in gqten's scale field.
    // Object metadata, allocator overhead and spare capacity are excluded.
    return tcapi::size(ctx, a) * sizeof(elem_t<TenT>);
}
template<class TenT> elem_t<TenT> get_elem(context_handle_t<TenT>& ctx, const TenT& a,
                                        const elem_coors_t<TenT>& coors) {
    detail::verbose::call diagnostic("get_elem", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "coors", coors);
    });
    detail::require_context(ctx);
    detail::check_coordinates(a, coors);
    return a.GetElem(coors);
}
}
