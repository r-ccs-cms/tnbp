#pragma once
#include "miscellaneous.h"

namespace tcapi {
namespace detail {
inline void check_inverse_status(int info) {
    if (info>0) throw std::domain_error("min-tcapi: singular matrix in inverse");
    if (info<0) throw std::runtime_error("min-tcapi: LAPACK inverse rejected an argument");
}
template<class E> void checked_matrix_inverse(int n, E* data) {
    using namespace gqten::hp_numeric;
    std::vector<int> pivots(static_cast<std::size_t>(n));
    int info=0;
    if constexpr(std::is_same_v<E,float>) sgetrf_(&n,&n,data,&n,pivots.data(),&info);
    else if constexpr(std::is_same_v<E,double>) dgetrf_(&n,&n,data,&n,pivots.data(),&info);
    else if constexpr(std::is_same_v<E,std::complex<float>>) cgetrf_(&n,&n,data,&n,pivots.data(),&info);
    else zgetrf_(&n,&n,data,&n,pivots.data(),&info);
    check_inverse_status(info);
    auto getri=[&](E* work,int length) {
        if constexpr(std::is_same_v<E,float>) sgetri_(&n,data,&n,pivots.data(),work,&length,&info);
        else if constexpr(std::is_same_v<E,double>) dgetri_(&n,data,&n,pivots.data(),work,&length,&info);
        else if constexpr(std::is_same_v<E,std::complex<float>>) cgetri_(&n,data,&n,pivots.data(),work,&length,&info);
        else zgetri_(&n,data,&n,pivots.data(),work,&length,&info);
        check_inverse_status(info);
    };
    E query{};getri(&query,-1);
    const long double requested=std::ceil(static_cast<long double>(std::real(query)));
    if (!std::isfinite(requested) || requested<n || requested>std::numeric_limits<int>::max())
        throw std::overflow_error("min-tcapi: invalid LAPACK inverse workspace size");
    const auto length=static_cast<int>(requested);
    std::vector<E> work(static_cast<std::size_t>(length));
    getri(work.data(),length);
}
} // namespace detail

template<class TenT>
void inverse(context_handle_t<TenT>& ctx, const TenT& in, order_t<TenT> rows, TenT& out) {
    detail::verbose::call diagnostic("inverse", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "in", in);
        detail::verbose::field(log, "rows", rows);
    });
    detail::require_context(ctx);
    if (rows<1 || static_cast<std::size_t>(rows)>=in.Rank())
        throw std::invalid_argument("min-tcapi: inverse requires 1 <= row bonds < rank");
    const auto count=detail::checked_size<TenT>(in.Shape());
    if (count>static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("min-tcapi: inverse exceeds backend integer size limit");
    std::size_t n=1;
    for (order_t<TenT> i=0;i<rows;++i) n*=in.Shape()[i];
    if (n!=count/n) throw std::invalid_argument("min-tcapi: inverse requires square matricization");
    std::vector<elem_t<TenT>> data(count);
    for (std::size_t i=0;i<count;++i) {
        data[i]=detail::numeric_value(in,i);
        if (!std::isfinite(std::real(data[i])) || !std::isfinite(std::imag(data[i])))
            throw std::domain_error("min-tcapi: inverse requires finite input");
    }
    detail::checked_matrix_inverse(static_cast<int>(n),data.data());
    auto result=detail::construct<TenT>(in.Shape(),[&](std::size_t i){return data[i];});
    detail::replace(out,std::move(result));
}
template<class TenT>
void inverse(context_handle_t<TenT>& ctx, TenT& inout, order_t<TenT> rows) {
    detail::verbose::call diagnostic("inverse", [&](std::ostream& log) {
        detail::verbose::field(log, "dtype", detail::verbose::scalar_name<elem_t<TenT>>());
        detail::verbose::field(log, "inout", inout);
        detail::verbose::field(log, "rows", rows);
    });
    tcapi::inverse(ctx,static_cast<const TenT&>(inout),rows,inout);
}
} // namespace tcapi
