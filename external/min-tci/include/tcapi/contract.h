#pragma once
#include "miscellaneous.h"
#include "gqten/tensor/out_of_place_ops/contract.h"
#include <algorithm>
#include <string_view>

namespace tcapi {
template<class TenT>
void contract(context_handle_t<TenT>& ctx,
              const TenT& a, const List<bond_label_t<TenT>>& la,
              const TenT& b, const List<bond_label_t<TenT>>& lb,
              TenT& out, const List<bond_label_t<TenT>>& lc) {
    detail::verbose::call diagnostic("contract", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "la", la);
        detail::verbose::field(log, "b", b);
        detail::verbose::field(log, "lb", lb);
        detail::verbose::field(log, "lc", lc);
    });
    detail::require_context(ctx);
    if (la.size()!=a.Rank() || lb.size()!=b.Rank())
        throw std::invalid_argument("min-tcapi: contract label count must match rank");
    for (const auto* labels : {&la,&lb,&lc}) {
        auto sorted=*labels;
        std::sort(sorted.begin(),sorted.end());
        if (std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end())
            throw std::invalid_argument("min-tcapi: repeated contract labels are unsupported");
        if (labels->size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::overflow_error("min-tcapi: too many contract labels");
    }
    const auto na=detail::checked_size<TenT>(a.Shape());
    const auto nb=detail::checked_size<TenT>(b.Shape());
    List<bond_label_t<TenT>> shared, free_labels;
    for (std::size_t i=0;i<la.size();++i) {
        auto j=std::find(lb.begin(),lb.end(),la[i]);
        if (j==lb.end()) free_labels.push_back(la[i]);
        else {
            if (a.Shape()[i]!=b.Shape()[j-lb.begin()])
                throw std::invalid_argument("min-tcapi: contracted dimensions differ");
            if (std::find(lc.begin(),lc.end(),la[i])!=lc.end())
                throw std::invalid_argument("min-tcapi: retaining shared contract labels is unsupported");
            shared.push_back(la[i]);
        }
    }
    for (auto label:lb)
        if (std::find(la.begin(),la.end(),label)==la.end()) free_labels.push_back(label);
    auto sorted_output=lc;
    std::sort(free_labels.begin(),free_labels.end());
    std::sort(sorted_output.begin(),sorted_output.end());
    if (free_labels!=sorted_output)
        throw std::invalid_argument("min-tcapi: output must contain every free label exactly once");
    shape_t<TenT> output_shape;
    for (auto label:lc) {
        auto i=std::find(la.begin(),la.end(),label);
        output_shape.push_back(i!=la.end() ? a.Shape()[i-la.begin()]
            : b.Shape()[std::find(lb.begin(),lb.end(),label)-lb.begin()]);
    }
    const auto nc=detail::checked_size<TenT>(output_shape);
    // gqten multiplies m*n, m*k and k*n as signed int, including allocations.
    const auto limit=static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (na>limit || nb>limit || nc>limit)
        throw std::overflow_error("min-tcapi: contract exceeds gqten integer size limit");
    TenT result;
    if (la.empty() || lb.empty()) {
        // Native Contract requires dense input storage, which scalars do not own.
        const auto& dense=la.empty()?b:a;
        const auto& labels=la.empty()?lb:la;
        const auto scalar=la.empty()?a.GetScale():b.GetScale();
        result=detail::construct<TenT>(output_shape,[&](std::size_t index) {
            auto coors=detail::coordinates<TenT>(index,output_shape);
            std::size_t source=0, stride=1;
            for (std::size_t i=0;i<labels.size();++i) {
                source+=coors[std::find(lc.begin(),lc.end(),labels[i])-lc.begin()]*stride;
                stride*=static_cast<std::size_t>(dense.Shape()[i]);
            }
            return scalar*detail::numeric_value(dense,source);
        });
    } else {
        // As in the legacy adapter: output positions are nonnegative;
        // contracted labels are mapped to -1, -2, ... independently of user IDs.
        std::sort(shared.begin(),shared.end());
        auto native_labels=[&](const auto& labels) {
            std::vector<int> mapped;
            for (auto label:labels) {
                auto i=std::find(shared.begin(),shared.end(),label);
                mapped.push_back(i!=shared.end() ? -1-static_cast<int>(i-shared.begin())
                    : static_cast<int>(std::find(lc.begin(),lc.end(),label)-lc.begin()));
            }
            return mapped;
        };
        gqten::Contract(&a,&b,native_labels(la),native_labels(lb),&result);
    }
    // Finish reading both inputs before replacing output, including out==a==b.
    detail::replace(out,std::move(result));
}

template<class TenT>
void contract(context_handle_t<TenT>& ctx, const TenT& a, std::string_view la,
              const TenT& b, std::string_view lb, TenT& out, std::string_view lc) {
    detail::verbose::call diagnostic("contract", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "a", a);
        detail::verbose::field(log, "la", la);
        detail::verbose::field(log, "b", b);
        detail::verbose::field(log, "lb", lb);
        detail::verbose::field(log, "lc", lc);
    });
    auto labels=[](std::string_view text) {
        List<bond_label_t<TenT>> result;
        for (unsigned char ch:text) result.push_back(ch);
        return result;
    };
    tcapi::contract(ctx,a,labels(la),b,labels(lb),out,labels(lc));
}
} // namespace tcapi
