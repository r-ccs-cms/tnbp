#pragma once
#include "read_only_queries.h"
#include "construction_and_destruction.h"
#include <cmath>
namespace tcapi {
inline void create_context(gqten_handle& ctx) {
    detail::verbose::call diagnostic("create_context", [](std::ostream&) {});
    ctx.active = true;
}
inline void destroy_context(gqten_handle& ctx) {
    detail::verbose::call diagnostic("destroy_context", [](std::ostream&) {});
    ctx.active = false;
}
template<class TenT, class RandomIt, class Func>
void to_range(context_handle_t<TenT>& ctx, const TenT& a, RandomIt first, Func&& coors2idx) {
    detail::verbose::call diagnostic("to_range", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    const auto count = tcapi::size(ctx, a);
    for (std::size_t i = 0; i < count; ++i) {
        const auto coors = detail::coordinates<TenT>(i, a.Shape());
        const auto index = detail::verbose::invoke_callback(coors2idx, coors);
        if (index < 0 || static_cast<std::size_t>(index) >= count)
            throw std::out_of_range("min-tcapi: output range index outside tensor size");
        *(first + index) = a.GetElem(coors);
    }
}

namespace detail {
template<class TenT> elem_t<TenT> numeric_value(const TenT& a, std::size_t i) {
    if (a.Rank() == 0) return a.GetScale();
    // Avoid multiplying complex nonfinite values by (1,0), which can corrupt
    // the other component. Non-unit lazy scales follow gqten arithmetic.
    if (a.GetScale() == elem_t<TenT>(1)) return a.GetRaw()[i];
    return a.GetElem(i);
}
template<class R, class S> R convert_component(S value) {
    // Define narrowing overflow explicitly instead of an out-of-range cast.
    if constexpr (sizeof(R) < sizeof(S)) {
        if (value > static_cast<S>(std::numeric_limits<R>::max())) return std::numeric_limits<R>::infinity();
        if (value < -static_cast<S>(std::numeric_limits<R>::max())) return -std::numeric_limits<R>::infinity();
    }
    return static_cast<R>(value);
}
}
template<class TenT>
bool close(context_handle_t<TenT>& ctx, const TenT& a, const TenT& b, const real_t<TenT> epsilon) {
    detail::verbose::call diagnostic("close", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "b", b);
        detail::verbose::field(log, "epsilon", epsilon);
    });
    detail::require_context(ctx);
    const auto count = detail::checked_size<TenT>(a.Shape());
    detail::checked_size<TenT>(b.Shape());
    if (!(epsilon >= real_t<TenT>(0))) throw std::invalid_argument("min-tcapi: tolerance must be nonnegative and not NaN");
    if (a.Shape() != b.Shape()) return false;
    for (std::size_t i=0; i<count; ++i) {
        const auto x = detail::numeric_value(a,i);
        const auto y = detail::numeric_value(b,i);
        real_t<TenT> difference;
        if constexpr (std::is_same_v<elem_t<TenT>,real_t<TenT>>) {
            if (!std::isfinite(x) || !std::isfinite(y)) return false;
            difference = std::abs(x-y);
        } else {
            if (!std::isfinite(x.real()) || !std::isfinite(x.imag()) ||
                !std::isfinite(y.real()) || !std::isfinite(y.imag())) return false;
            difference = std::hypot(x.real()-y.real(),x.imag()-y.imag());
        }
        if (!(difference <= epsilon)) return false;
    }
    return true;
}
template<class Ten1T, class Ten2T>
void convert(context_handle_t<Ten1T>& ctx1, const Ten1T& in,
             context_handle_t<Ten2T>& ctx2, Ten2T& out) {
    detail::verbose::call diagnostic("convert", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<Ten2T>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx1);
    detail::require_context(ctx2);
    auto result = [&]() -> Ten2T {
        if constexpr (std::is_same_v<Ten1T,Ten2T>) return tcapi::copy(ctx2,in);
        else {
            auto converted = detail::construct<Ten2T>(in.Shape(),[&](std::size_t i) {
                const auto e = detail::numeric_value(in,i);
                using R = real_t<Ten2T>;
                const auto r = detail::convert_component<R>(std::real(e));
                if constexpr (std::is_same_v<elem_t<Ten2T>,R>) return r;
                else return elem_t<Ten2T>(r,detail::convert_component<R>(std::imag(e)));
            });
            if (in.Labeled()) converted.SetLabels(in.GetLabels());
            return converted;
        }
    }();
    detail::replace(out,std::move(result));
}
}

namespace tcapi {
template<class TenT> std::string version() {
    detail::verbose::call diagnostic("version", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
    });
    return "1.0"; // User-selected specification baseline; limitations are documented.
}
template<class TenT>
void show(context_handle_t<TenT>& ctx, const TenT& a) {
    detail::verbose::call diagnostic("show", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(a.Shape());
    a.FormattedPrint();
}
}
