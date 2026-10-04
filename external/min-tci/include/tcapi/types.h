#pragma once
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gqten { template<class T> class tensor; }
namespace tcapi {
template<class T> using List = std::vector<T>;
template<class T> using CRef = std::reference_wrapper<const T>;
template<class T, class U> using Pair = std::pair<T, U>;
template<class T, class U> using Map = std::unordered_map<T, U>;
struct gqten_handle { bool active = false; };
namespace detail {
template<class T> struct scalar_types;
template<> struct scalar_types<float> { using real = float; using complex = std::complex<float>; };
template<> struct scalar_types<double> { using real = double; using complex = std::complex<double>; };
template<> struct scalar_types<std::complex<float>> : scalar_types<float> {};
template<> struct scalar_types<std::complex<double>> : scalar_types<double> {};
}
template<class TenT> struct tensor_traits;
template<class T> struct tensor_traits<gqten::tensor<T>> {
    using ten_t = gqten::tensor<T>;
    using order_t = std::int32_t;
    using bond_dim_t = std::int32_t;
    using bond_idx_t = std::int32_t;
    using bond_label_t = std::int32_t;
    using ten_size_t = std::size_t;
    using elem_t = T;
    using elem_coor_t = std::int32_t;
    using shape_t = List<bond_dim_t>;
    using elem_coors_t = List<elem_coor_t>;
    using real_t = typename detail::scalar_types<T>::real;
    using real_ten_t = gqten::tensor<real_t>;
    using cplx_t = typename detail::scalar_types<T>::complex;
    using cplx_ten_t = gqten::tensor<cplx_t>;
    using context_handle_t = gqten_handle;
};
#define TCAPI_ASSOCIATED_TYPE(name) template<class T> using name = typename tensor_traits<T>::name;
TCAPI_ASSOCIATED_TYPE(ten_t)
TCAPI_ASSOCIATED_TYPE(order_t)
TCAPI_ASSOCIATED_TYPE(bond_dim_t)
TCAPI_ASSOCIATED_TYPE(bond_idx_t)
TCAPI_ASSOCIATED_TYPE(bond_label_t)
TCAPI_ASSOCIATED_TYPE(ten_size_t)
TCAPI_ASSOCIATED_TYPE(elem_t)
TCAPI_ASSOCIATED_TYPE(elem_coor_t)
TCAPI_ASSOCIATED_TYPE(shape_t)
TCAPI_ASSOCIATED_TYPE(elem_coors_t)
TCAPI_ASSOCIATED_TYPE(real_t)
TCAPI_ASSOCIATED_TYPE(real_ten_t)
TCAPI_ASSOCIATED_TYPE(cplx_t)
TCAPI_ASSOCIATED_TYPE(cplx_ten_t)
TCAPI_ASSOCIATED_TYPE(context_handle_t)
#undef TCAPI_ASSOCIATED_TYPE
}
