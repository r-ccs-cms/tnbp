#pragma once
#include "construction_and_destruction.h"
#include <algorithm>

namespace tcapi {
namespace detail {
inline void check_axis(std::int32_t axis, std::size_t rank) {
    if (axis < 0 || static_cast<std::size_t>(axis) >= rank)
        throw std::out_of_range("min-tcapi: bond index outside tensor");
}
template<class TenT>
TenT slice(const TenT& in, const List<Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>>& ranges) {
    checked_size<TenT>(in.Shape());
    if (ranges.size() != in.Rank()) throw std::invalid_argument("min-tcapi: slice rank mismatch");
    auto shape = in.Shape();
    for (std::size_t k=0; k<shape.size(); ++k) {
        const auto [first,last] = ranges[k];
        if (first < 0 || last < first || last > shape[k])
            throw std::out_of_range("min-tcapi: invalid half-open slice range");
        shape[k] = last-first; // checked_size rejects empty slices before allocation.
    }
    auto result = construct<TenT>(shape, [&](std::size_t i) {
        auto c = coordinates<TenT>(i,shape);
        for (std::size_t k=0; k<c.size(); ++k) c[k] += ranges[k].first;
        return in.GetElem(c);
    });
    if (in.Labeled()) result.SetLabels(in.GetLabels());
    return result;
}
}
template<class TenT>
void expand(context_handle_t<TenT>& ctx, const TenT& in,
            const Map<bond_idx_t<TenT>,bond_dim_t<TenT>>& increments, TenT& out) {
    detail::verbose::call diagnostic("expand", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "increments", increments);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(in.Shape());
    auto shape = in.Shape();
    for (const auto& [axis,increment] : increments) {
        detail::check_axis(axis,shape.size());
        if (increment < 0) throw std::invalid_argument("min-tcapi: negative expansion");
        if (increment > std::numeric_limits<bond_dim_t<TenT>>::max()-shape[axis])
            throw std::overflow_error("min-tcapi: expanded dimension overflow");
        shape[axis] += increment;
    }
    auto result = detail::construct<TenT>(shape, [&](std::size_t i) {
        const auto c = detail::coordinates<TenT>(i,shape);
        for (std::size_t k=0; k<c.size(); ++k)
            if (c[k] >= in.Shape()[k]) return elem_t<TenT>{};
        return in.GetElem(c);
    });
    if (in.Labeled()) result.SetLabels(in.GetLabels());
    detail::replace(out,std::move(result));
}
template<class TenT>
void expand(context_handle_t<TenT>& ctx, TenT& inout,
            const Map<bond_idx_t<TenT>,bond_dim_t<TenT>>& increments) {
    detail::verbose::call diagnostic("expand", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "increments", increments);
    });
    tcapi::expand(ctx,static_cast<const TenT&>(inout),increments,inout);
}
template<class TenT>
void shrink(context_handle_t<TenT>& ctx, const TenT& in,
            const Map<bond_idx_t<TenT>,Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>>& ranges, TenT& out) {
    detail::verbose::call diagnostic("shrink", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "ranges", ranges);
    });
    detail::require_context(ctx);
    List<Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>> all;
    for (auto dim : in.Shape()) all.emplace_back(0,dim);
    for (const auto& [axis,range] : ranges) {
        detail::check_axis(axis,all.size());
        all[axis] = range;
    }
    auto result = detail::slice(in,all);
    detail::replace(out,std::move(result));
}
template<class TenT>
void shrink(context_handle_t<TenT>& ctx, TenT& inout,
            const Map<bond_idx_t<TenT>,Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>>& ranges) {
    detail::verbose::call diagnostic("shrink", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "ranges", ranges);
    });
    tcapi::shrink(ctx,static_cast<const TenT&>(inout),ranges,inout);
}
template<class TenT>
void extract_sub(context_handle_t<TenT>& ctx, const TenT& in,
                 const List<Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>>& ranges, TenT& out) {
    detail::verbose::call diagnostic("extract_sub", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "ranges", ranges);
    });
    detail::require_context(ctx);
    auto result = detail::slice(in,ranges);
    detail::replace(out,std::move(result));
}
template<class TenT>
void extract_sub(context_handle_t<TenT>& ctx, TenT& inout,
                 const List<Pair<elem_coor_t<TenT>,elem_coor_t<TenT>>>& ranges) {
    detail::verbose::call diagnostic("extract_sub", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "ranges", ranges);
    });
    tcapi::extract_sub(ctx,static_cast<const TenT&>(inout),ranges,inout);
}
template<class TenT>
void replace_sub(context_handle_t<TenT>& ctx, const TenT& in, const TenT& sub,
                 const elem_coors_t<TenT>& begin, TenT& out) {
    detail::verbose::call diagnostic("replace_sub", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "sub", sub);
        detail::verbose::field(log, "begin", begin);
    });
    detail::require_context(ctx);
    detail::checked_size<TenT>(in.Shape());
    detail::checked_size<TenT>(sub.Shape());
    if (in.Rank() != sub.Rank() || begin.size() != in.Rank())
        throw std::invalid_argument("min-tcapi: replacement rank mismatch");
    for (std::size_t k=0; k<begin.size(); ++k)
        if (begin[k] < 0 || begin[k] > in.Shape()[k] || sub.Shape()[k] > in.Shape()[k]-begin[k])
            throw std::out_of_range("min-tcapi: replacement outside tensor");
    auto result = detail::construct<TenT>(in.Shape(), [&](std::size_t i) {
        auto c = detail::coordinates<TenT>(i,in.Shape());
        auto local = c;
        bool inside = true;
        for (std::size_t k=0; k<c.size(); ++k) {
            local[k] -= begin[k];
            if (local[k] < 0 || local[k] >= sub.Shape()[k]) inside = false;
        }
        return inside ? sub.GetElem(local) : in.GetElem(c);
    });
    if (in.Labeled()) result.SetLabels(in.GetLabels());
    detail::replace(out,std::move(result));
}
template<class TenT>
void replace_sub(context_handle_t<TenT>& ctx, TenT& inout, const TenT& sub,
                 const elem_coors_t<TenT>& begin) {
    detail::verbose::call diagnostic("replace_sub", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "sub", sub);
        detail::verbose::field(log, "begin", begin);
    });
    tcapi::replace_sub(ctx,static_cast<const TenT&>(inout),sub,begin,inout);
}
template<class TenT>
TenT concatenate(context_handle_t<TenT>& ctx, const List<CRef<TenT>>& ins, const bond_idx_t<TenT> axis) {
    detail::verbose::call diagnostic("concatenate", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "ins", ins);
        detail::verbose::field(log, "axis", axis);
    });
    detail::require_context(ctx);
    if (ins.empty()) throw std::invalid_argument("min-tcapi: empty concatenation input");
    auto shape = ins.front().get().Shape();
    detail::check_axis(axis,shape.size());
    shape[axis] = 0;
    std::vector<bond_dim_t<TenT>> ends;
    ends.reserve(ins.size());
    for (const auto& ref : ins) {
        const auto& in = ref.get();
        detail::checked_size<TenT>(in.Shape());
        if (in.Rank() != shape.size()) throw std::invalid_argument("min-tcapi: concatenation rank mismatch");
        for (std::size_t k=0; k<shape.size(); ++k)
            if (k != static_cast<std::size_t>(axis) && in.Shape()[k] != shape[k])
                throw std::invalid_argument("min-tcapi: concatenation shape mismatch");
        if (in.Shape()[axis] > std::numeric_limits<bond_dim_t<TenT>>::max()-shape[axis])
            throw std::overflow_error("min-tcapi: concatenated dimension overflow");
        shape[axis] += in.Shape()[axis];
        ends.push_back(shape[axis]);
    }
    return detail::construct<TenT>(shape,[&](std::size_t i) {
        auto c = detail::coordinates<TenT>(i,shape);
        const auto index = static_cast<std::size_t>(std::upper_bound(ends.begin(),ends.end(),c[axis])-ends.begin());
        if (index) c[axis] -= ends[index-1];
        return ins[index].get().GetElem(c);
    });
}
template<class TenT>
TenT stack(context_handle_t<TenT>& ctx, const List<CRef<TenT>>& ins, const bond_idx_t<TenT> axis) {
    detail::verbose::call diagnostic("stack", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "ins", ins);
        detail::verbose::field(log, "axis", axis);
    });
    detail::require_context(ctx);
    if (ins.empty()) throw std::invalid_argument("min-tcapi: empty stack input");
    auto shape = ins.front().get().Shape();
    if (axis < 0 || static_cast<std::size_t>(axis) > shape.size())
        throw std::out_of_range("min-tcapi: stack axis outside insertion range");
    if (ins.size() > static_cast<std::size_t>(std::numeric_limits<bond_dim_t<TenT>>::max()))
        throw std::overflow_error("min-tcapi: stack dimension overflow");
    for (const auto& ref : ins) {
        detail::checked_size<TenT>(ref.get().Shape());
        if (ref.get().Shape() != shape) throw std::invalid_argument("min-tcapi: stack shape mismatch");
    }
    shape.insert(shape.begin()+axis,static_cast<bond_dim_t<TenT>>(ins.size()));
    return detail::construct<TenT>(shape,[&](std::size_t i) {
        auto c = detail::coordinates<TenT>(i,shape);
        const auto index = static_cast<std::size_t>(c[axis]);
        c.erase(c.begin()+axis);
        return ins[index].get().GetElem(c);
    });
}
}
