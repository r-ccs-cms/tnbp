#pragma once
#include "construction_and_destruction.h"
namespace tcapi {
template<class TenT>
void set_elem(context_handle_t<TenT>& ctx, TenT& a, const elem_coors_t<TenT>& coors,
              const elem_t<TenT> value) {
    detail::verbose::call diagnostic("set_elem", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "coors", coors);
        detail::verbose::field(log, "value", value);
    });
    detail::require_context(ctx);
    detail::check_coordinates(a, coors);
    // gqten SetElem silently ignores dense writes when its lazy scale is zero.
    if (a.Rank() != 0 && a.GetScale() == elem_t<TenT>{}) {
        auto replacement = tcapi::zeros<TenT>(ctx, a.Shape());
        replacement.SetElem(coors, value);
        a = std::move(replacement);
    } else {
        a.SetElem(coors, value);
    }
}

template<class TenT>
void reshape(context_handle_t<TenT>& ctx, TenT& inout, const shape_t<TenT>& new_shape) {
    detail::verbose::call diagnostic("reshape", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "new_shape", new_shape);
    });
    detail::require_context(ctx);
    if (detail::checked_size<TenT>(new_shape) != detail::checked_size<TenT>(inout.Shape()))
        throw std::invalid_argument("min-tcapi: reshape changes the logical size");
    if (inout.Rank() == 0 || new_shape.empty()) {
        const auto value = inout.Rank() == 0 ? inout.GetScale() : inout.GetElem(std::size_t{0});
        auto result = detail::construct<TenT>(new_shape, [value](std::size_t) { return value; });
        detail::replace(inout, std::move(result));
    } else {
        inout.Reshape(new_shape);
    }
}
template<class TenT>
void reshape(context_handle_t<TenT>& ctx, const TenT& in, const shape_t<TenT>& new_shape, TenT& out) {
    detail::verbose::call diagnostic("reshape", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "new_shape", new_shape);
    });
    detail::require_context(ctx);
    if (detail::checked_size<TenT>(new_shape) != detail::checked_size<TenT>(in.Shape()))
        throw std::invalid_argument("min-tcapi: reshape changes the logical size");
    auto result = tcapi::copy(ctx, in);
    tcapi::reshape(ctx, result, new_shape);
    detail::replace(out, std::move(result));
}

template<class TenT>
void transpose(context_handle_t<TenT>& ctx, TenT& inout,
               const List<bond_idx_t<TenT>>& new_order) {
    detail::verbose::call diagnostic("transpose", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "new_order", new_order);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(inout.Shape());
    if (new_order.size() != inout.Rank())
        throw std::invalid_argument("min-tcapi: permutation rank mismatch");
    std::vector<bool> seen(new_order.size(), false);
    for (const auto axis : new_order) {
        if (axis < 0 || static_cast<std::size_t>(axis) >= new_order.size())
            throw std::out_of_range("min-tcapi: permutation axis outside tensor");
        if (seen[axis]) throw std::invalid_argument("min-tcapi: duplicate permutation axis");
        seen[axis] = true;
    }
    inout.Transpose(new_order);
}
template<class TenT>
void transpose(context_handle_t<TenT>& ctx, const TenT& in,
               const List<bond_idx_t<TenT>>& new_order, TenT& out) {
    detail::verbose::call diagnostic("transpose", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "new_order", new_order);
    });
    auto result = tcapi::copy(ctx, in);
    tcapi::transpose(ctx, result, new_order);
    detail::replace(out, std::move(result));
}

template<class TenT>
void cplx_conj(context_handle_t<TenT>& ctx, TenT& inout) {
    detail::verbose::call diagnostic("cplx_conj", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
    });
    detail::require_context(ctx);
    if constexpr (!std::is_same_v<elem_t<TenT>, real_t<TenT>>) {
        if (inout.Rank() == 0) inout.SetScale(std::conj(inout.GetScale()));
        else inout.Conjugate();
    }
}
template<class TenT>
void cplx_conj(context_handle_t<TenT>& ctx, const TenT& in, TenT& out) {
    detail::verbose::call diagnostic("cplx_conj", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    auto result = tcapi::copy(ctx, in);
    tcapi::cplx_conj(ctx, result);
    detail::replace(out, std::move(result));
}

namespace detail {
template<class TenT> elem_t<TenT> logical_element(const TenT& in, std::size_t i) {
    return in.Rank() == 0 ? in.GetScale() : in.GetElem(i);
}
template<class OutT, class InT, class Func>
OutT map_elements(const InT& in, Func&& f) {
    auto result = construct<OutT>(in.Shape(), [&](std::size_t i) {
        return std::invoke(f, logical_element(in, i));
    });
    if (in.Labeled()) result.SetLabels(in.GetLabels());
    return result;
}
}
template<class TenT>
cplx_ten_t<TenT> to_cplx(context_handle_t<TenT>& ctx, const TenT& in) {
    detail::verbose::call diagnostic("to_cplx", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx);
    if constexpr (std::is_same_v<elem_t<TenT>, cplx_t<TenT>>) return tcapi::copy(ctx, in);
    else return detail::map_elements<cplx_ten_t<TenT>>(in, [](elem_t<TenT> e) {
        return cplx_t<TenT>(e, 0);
    });
}
template<class TenT>
real_ten_t<TenT> real(context_handle_t<TenT>& ctx, const TenT& in) {
    detail::verbose::call diagnostic("real", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx);
    if constexpr (std::is_same_v<elem_t<TenT>, real_t<TenT>>) return tcapi::copy(ctx, in);
    else return detail::map_elements<real_ten_t<TenT>>(in, [](elem_t<TenT> e) { return e.real(); });
}
template<class TenT>
real_ten_t<TenT> imag(context_handle_t<TenT>& ctx, const TenT& in) {
    detail::verbose::call diagnostic("imag", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx);
    if constexpr (std::is_same_v<elem_t<TenT>, real_t<TenT>>) {
        auto result = tcapi::zeros<real_ten_t<TenT>>(ctx, in.Shape());
        if (in.Labeled()) result.SetLabels(in.GetLabels());
        return result;
    } else return detail::map_elements<real_ten_t<TenT>>(in, [](elem_t<TenT> e) { return e.imag(); });
}

template<class TenT, class Func>
void for_each(context_handle_t<TenT>& ctx, TenT& inout, Func&& f) {
    detail::verbose::call diagnostic("for_each", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
    });
    detail::require_context(ctx);
    static_assert(std::is_invocable_v<Func&, elem_t<TenT>&>, "callback must accept a mutable element");
    auto result = detail::map_elements<TenT>(inout, [&](elem_t<TenT> e) {
        detail::verbose::invoke_callback(f, e);
        return e;
    });
    detail::replace(inout, std::move(result));
}
template<class TenT, class Func>
void for_each(context_handle_t<TenT>& ctx, const TenT& in, Func&& f) {
    detail::verbose::call diagnostic("for_each", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx);
    static_assert(std::is_invocable_v<Func&, const elem_t<TenT>&>, "callback must accept a const element");
    const auto count = detail::checked_size<TenT>(in.Shape());
    for (std::size_t i = 0; i < count; ++i) {
        const auto e = detail::logical_element(in, i);
        detail::verbose::invoke_callback(f, e);
    }
}
template<class TenT, class Func>
void for_each_with_coors(context_handle_t<TenT>& ctx, TenT& inout, Func&& f) {
    detail::verbose::call diagnostic("for_each_with_coors", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
    });
    detail::require_context(ctx);
    static_assert(std::is_invocable_v<Func&, elem_t<TenT>&, const elem_coors_t<TenT>&>,
                  "callback must accept a mutable element and const coordinates");
    auto result = detail::construct<TenT>(inout.Shape(), [&](std::size_t i) {
        const auto coors = detail::coordinates<TenT>(i, inout.Shape());
        auto e = detail::logical_element(inout, i);
        detail::verbose::invoke_callback(f, e, coors);
        return e;
    });
    if (inout.Labeled()) result.SetLabels(inout.GetLabels());
    detail::replace(inout, std::move(result));
}
template<class TenT, class Func>
void for_each_with_coors(context_handle_t<TenT>& ctx, const TenT& in, Func&& f) {
    detail::verbose::call diagnostic("for_each_with_coors", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx);
    static_assert(std::is_invocable_v<Func&, const elem_t<TenT>&, const elem_coors_t<TenT>&>,
                  "callback must accept a const element and const coordinates");
    const auto count = detail::checked_size<TenT>(in.Shape());
    for (std::size_t i = 0; i < count; ++i) {
        const auto coors = detail::coordinates<TenT>(i, in.Shape());
        const auto e = detail::logical_element(in, i);
        detail::verbose::invoke_callback(f, e, coors);
    }
}
}

#include "tensor_regions.h"
