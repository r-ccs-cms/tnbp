#pragma once
#include "miscellaneous.h"
#include <algorithm>
#include "gqten/tensor/out_of_place_ops/linear_combine.h"

namespace tcapi {
template<class TenT>
void diag(context_handle_t<TenT>& ctx, const TenT& in, TenT& out) {
    detail::verbose::call diagnostic("diag", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(in.Shape());
    if (in.Rank() != 1 && in.Rank() != 2)
        throw std::invalid_argument("min-tcapi: diag requires a vector or matrix");
    const auto n = in.Rank() == 1 ? in.Shape()[0] : std::min(in.Shape()[0],in.Shape()[1]);
    const shape_t<TenT> shape = in.Rank() == 1 ? shape_t<TenT>{n,n} : shape_t<TenT>{n};
    auto result = detail::construct<TenT>(shape,[&](std::size_t i) {
        if (in.Rank() == 1) {
            const auto dim = static_cast<std::size_t>(n);
            return i % dim == i / dim ? detail::numeric_value(in,i % dim) : elem_t<TenT>{};
        }
        return detail::numeric_value(in,i*(static_cast<std::size_t>(in.Shape()[0])+1));
    });
    detail::replace(out,std::move(result));
}
template<class TenT>
void diag(context_handle_t<TenT>& ctx, TenT& inout) {
    detail::verbose::call diagnostic("diag", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
    });
    tcapi::diag(ctx,static_cast<const TenT&>(inout),inout);
}
template<class TenT>
real_t<TenT> norm(context_handle_t<TenT>& ctx, const TenT& a) {
    detail::verbose::call diagnostic("norm", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
    });
    detail::require_context(ctx);
    const auto count = detail::checked_size<TenT>(a.Shape());
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("min-tcapi: gqten BLAS vector length exceeds int");
    return a.CalcNorm();
}
template<class TenT>
void scale(context_handle_t<TenT>& ctx, TenT& inout, const elem_t<TenT> s) {
    detail::verbose::call diagnostic("scale", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "s", s);
    });
    detail::require_context(ctx);
    inout.Scale(s);
}
template<class TenT>
void scale(context_handle_t<TenT>& ctx, const TenT& in, const elem_t<TenT> s, TenT& out) {
    detail::verbose::call diagnostic("scale", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "s", s);
    });
    if (std::addressof(in) == std::addressof(out)) { tcapi::scale(ctx,out,s); return; }
    auto result = tcapi::copy(ctx,in);
    tcapi::scale(ctx,result,s);
    detail::replace(out,std::move(result));
}
namespace detail {
template<class TenT> void normalize_by(TenT& a, real_t<TenT> n) {
    const auto factor = a.GetScale()/n;
    if (std::isfinite(std::real(factor)) && std::isfinite(std::imag(factor)) && factor != elem_t<TenT>{}) {
        a.SetScale(factor);
    } else {
        // A reciprocal encoded only in the lazy scale can overflow/underflow
        // even when the normalized logical elements are representable.
        auto result = construct<TenT>(a.Shape(),[&](std::size_t i) { return numeric_value(a,i)/n; });
        if (a.Labeled()) result.SetLabels(a.GetLabels());
        replace(a,std::move(result));
    }
}
}
template<class TenT>
real_t<TenT> normalize(context_handle_t<TenT>& ctx, TenT& inout) {
    detail::verbose::call diagnostic("normalize", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
    });
    const auto n = tcapi::norm(ctx,inout);
    if (!(n > real_t<TenT>(0)) || !std::isfinite(n))
        throw std::domain_error("min-tcapi: normalization requires a finite positive norm");
    detail::normalize_by(inout,n);
    return n;
}
template<class TenT>
real_t<TenT> normalize(context_handle_t<TenT>& ctx, const TenT& in, TenT& out) {
    detail::verbose::call diagnostic("normalize", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
    });
    if (std::addressof(in) == std::addressof(out)) return tcapi::normalize(ctx,out);
    const auto n = tcapi::norm(ctx,in);
    if (!(n > real_t<TenT>(0)) || !std::isfinite(n))
        throw std::domain_error("min-tcapi: normalization requires a finite positive norm");
    auto result = tcapi::copy(ctx,in);
    detail::normalize_by(result,n);
    detail::replace(out,std::move(result));
    return n;
}

template<class TenT>
TenT linear_combine(context_handle_t<TenT>& ctx, const List<CRef<TenT>>& ins,
                    const List<elem_t<TenT>>& coefs) {
    detail::verbose::call diagnostic("linear_combine", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "ins", ins);
        detail::verbose::field(log, "coefs", coefs);
    });
    detail::require_context(ctx);
    if (ins.empty() || ins.size() != coefs.size())
        throw std::invalid_argument("min-tcapi: nonempty inputs and matching coefficients required");
    const auto& shape = ins.front().get().Shape();
    const auto count = detail::checked_size<TenT>(shape);
    for (const auto& ref : ins)
        if (ref.get().Shape() != shape) throw std::invalid_argument("min-tcapi: linear combination shape mismatch");
    if (shape.empty()) {
        elem_t<TenT> sum{};
        for (std::size_t i=0; i<ins.size(); ++i) sum += coefs[i]*ins[i].get().GetScale();
        return tcapi::fill<TenT>(ctx,{},sum);
    }
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("min-tcapi: gqten BLAS vector length exceeds int");
    std::vector<TenT*> inputs;
    inputs.reserve(ins.size());
    // gqten's legacy signature takes mutable pointers, but LinearCombine only
    // reads these inputs. Output is separate, including for repeated references.
    for (const auto& ref : ins) inputs.push_back(const_cast<TenT*>(std::addressof(ref.get())));
    auto result = tcapi::allocate<TenT>(ctx,shape);
    gqten::LinearCombine(coefs,inputs,&result);
    return result;
}
template<class TenT>
TenT linear_combine(context_handle_t<TenT>& ctx, const List<CRef<TenT>>& ins) {
    detail::verbose::call diagnostic("linear_combine", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "ins", ins);
    });
    return tcapi::linear_combine(ctx,ins,List<elem_t<TenT>>(ins.size(),elem_t<TenT>(1)));
}
template<class TenT>
void trace(context_handle_t<TenT>& ctx, const TenT& in,
           const List<Pair<bond_idx_t<TenT>,bond_idx_t<TenT>>>& pairs, TenT& out) {
    detail::verbose::call diagnostic("trace", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "pairs", pairs);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(in.Shape());
    std::vector<bool> used(in.Rank(),false);
    shape_t<TenT> traced_shape, output_shape;
    std::vector<std::size_t> remaining;
    for (const auto& [a,b] : pairs) {
        if (a < 0 || b < 0 || static_cast<std::size_t>(a) >= in.Rank() || static_cast<std::size_t>(b) >= in.Rank())
            throw std::out_of_range("min-tcapi: trace axis outside tensor");
        if (a == b || used[a] || used[b]) throw std::invalid_argument("min-tcapi: trace pairs must use distinct axes");
        if (in.Shape()[a] != in.Shape()[b]) throw std::invalid_argument("min-tcapi: trace dimensions differ");
        used[a] = used[b] = true;
        traced_shape.push_back(in.Shape()[a]);
    }
    if (pairs.empty()) {
        auto result = tcapi::copy(ctx,in);
        detail::replace(out,std::move(result));
        return;
    }
    for (std::size_t i=0; i<in.Rank(); ++i) if (!used[i]) {
        remaining.push_back(i); output_shape.push_back(in.Shape()[i]);
    }
    const auto traced_count = detail::checked_size<TenT>(traced_shape);
    // Same coordinate-sum method as gqten::Trace, without including trace.h:
    // that header defines a non-inline free helper and breaks multi-TU linking.
    auto result = detail::construct<TenT>(output_shape,[&](std::size_t i) {
        const auto output_coors = detail::coordinates<TenT>(i,output_shape);
        elem_coors_t<TenT> c(in.Rank(),0);
        for (std::size_t k=0; k<remaining.size(); ++k) c[remaining[k]] = output_coors[k];
        elem_t<TenT> sum{};
        for (std::size_t j=0; j<traced_count; ++j) {
            const auto t = detail::coordinates<TenT>(j,traced_shape);
            for (std::size_t k=0; k<pairs.size(); ++k) c[pairs[k].first] = c[pairs[k].second] = t[k];
            sum += in.GetElem(c);
        }
        return sum;
    });
    detail::replace(out,std::move(result));
}
template<class TenT>
void trace(context_handle_t<TenT>& ctx, TenT& inout,
           const List<Pair<bond_idx_t<TenT>,bond_idx_t<TenT>>>& pairs) {
    detail::verbose::call diagnostic("trace", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "pairs", pairs);
    });
    tcapi::trace(ctx,static_cast<const TenT&>(inout),pairs,inout);
}
}

#include "contract.h"
#include "svd.h"
#include "qr.h"
#include "eigh.h"
#include "exp.h"
#include "inverse.h"
