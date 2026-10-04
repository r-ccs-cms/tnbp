#pragma once
#include "detail/core.h"
#include <cstring>
#include <iterator>
namespace tcapi {
template<class TenT> TenT allocate(context_handle_t<TenT>& ctx, const shape_t<TenT>& shape) {
    detail::verbose::call diagnostic("allocate", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "shape", shape);
    });
    detail::require_context(ctx);
    const auto count = detail::checked_size<TenT>(shape);
    if (shape.empty()) return TenT{}; // The initial scalar value is unspecified by this API.
    using E = elem_t<TenT>;
    std::unique_ptr<E, decltype(&std::free)> data(
        static_cast<E*>(std::malloc(count * sizeof(E))), &std::free);
    if (!data) throw std::bad_alloc();
    TenT a(shape, data.get());
    data.release();
    return a;
}
template<class TenT> TenT fill(context_handle_t<TenT>& ctx, const shape_t<TenT>& shape, const elem_t<TenT> value) {
    detail::verbose::call diagnostic("fill", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "shape", shape);
        detail::verbose::field(log, "value", value);
    });
    detail::require_context(ctx);
    return detail::construct<TenT>(shape, [value](std::size_t) { return value; });
}
template<class TenT> TenT zeros(context_handle_t<TenT>& ctx, const shape_t<TenT>& shape) {
    detail::verbose::call diagnostic("zeros", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "shape", shape);
    });
    return tcapi::fill<TenT>(ctx, shape, elem_t<TenT>{});
}
template<class TenT> TenT eye(context_handle_t<TenT>& ctx, const bond_dim_t<TenT> n) {
    detail::verbose::call diagnostic("eye", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "n", n);
    });
    detail::require_context(ctx);
    return detail::construct<TenT>({n,n}, [n](std::size_t i) {
        const auto d = static_cast<std::size_t>(n);
        return elem_t<TenT>(i % d == i / d ? 1 : 0);
    });
}
template<class TenT, class RandomIt, class Func>
TenT assign_from_range(context_handle_t<TenT>& ctx, const shape_t<TenT>& shape,
                      RandomIt first, Func&& coors2idx) {
    detail::verbose::call diagnostic("assign_from_range", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "shape", shape);
    });
    detail::require_context(ctx);
    return detail::construct<TenT>(shape, [&](std::size_t i) {
        const auto coors = detail::coordinates<TenT>(i, shape);
        const auto index = detail::verbose::invoke_callback(coors2idx, coors);
        if (index < 0) throw std::out_of_range("min-tcapi: negative input range index");
        // The input range length is not supplied; its upper bound is the caller's responsibility.
        return *(first + index);
    });
}
template<class TenT, class RandNumGen>
TenT random(context_handle_t<TenT>& ctx, const shape_t<TenT>& shape, RandNumGen& gen) {
    detail::verbose::call diagnostic("random", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "shape", shape);
    });
    detail::require_context(ctx);
    return detail::construct<TenT>(shape, [&](std::size_t) { return detail::verbose::invoke_callback(gen); });
}
template<class TenT> TenT copy(context_handle_t<TenT>& ctx, const TenT& orig) {
    detail::verbose::call diagnostic("copy", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "orig", orig);
    });
    detail::require_context(ctx);
    auto result = tcapi::allocate<TenT>(ctx, orig.Shape());
    if (orig.Rank() != 0) {
        // Copy the representation without reading indeterminate elements or
        // rounding a lazy scale into the payload. result owns mutable storage.
        std::memcpy(const_cast<elem_t<TenT>*>(result.GetRaw()), orig.GetRaw(),
                    detail::checked_size<TenT>(orig.Shape()) * sizeof(elem_t<TenT>));
    }
    result.SetScale(orig.GetScale());
    if (orig.Labeled()) result.SetLabels(orig.GetLabels());
    return result;
}
template<class TenT> TenT move(context_handle_t<TenT>& ctx, TenT& from) {
    detail::verbose::call diagnostic("move", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "from", from);
    });
    detail::require_context(ctx);
    return TenT(std::move(from));
}
template<class TenT> void clear(context_handle_t<TenT>& ctx, TenT& a) {
    detail::verbose::call diagnostic("clear", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    detail::require_context(ctx);
    detail::reset(a);
}
}
