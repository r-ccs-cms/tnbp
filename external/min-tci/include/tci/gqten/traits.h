#ifndef TCI_GQTEN_TRAITS_H
#define TCI_GQTEN_TRAITS_H

// Aligned to the original traits (traits_original.h) exactly:
// - explicit specializations for float, double, complex<float>, complex<double>
// - shape_t, elem_coor_t, elem_coors_t, bond_* types match the original
// - provides gqten_handle and is_gqten_handle_v
// - minimal C++17-only traits; no concepts/requires

#include <vector>
#include <complex>
#include <cstdint>
#include <type_traits>
#include <string>

namespace gqten { template <typename ElemT> class tensor; }

namespace tci {

  template<class T> struct is_gqten_tensor : std::false_type {};
  template<class E> struct is_gqten_tensor<gqten::tensor<E>> : std::true_type {};
  template<class T>
  inline constexpr bool is_gqten_tensor_v = is_gqten_tensor<T>::value;

  struct gqten_handle {
    std::string text = std::string("none");
  };
  
  template <class T> struct is_gqten_handle : std::false_type {};
  template <> struct is_gqten_handle<gqten_handle> : std::true_type {};
  template <class T>
  inline constexpr bool is_gqten_handle_v =
    is_gqten_handle<std::remove_cv_t<std::remove_reference_t<T>>>::value;

  // traits for gqten::tensor<float>
  template <>
  struct tensor_traits<gqten::tensor<float>> {
    using ten_t = gqten::tensor<float>;
    using order_t = int32_t;
    using shape_t = std::vector<int32_t>;
    using bond_dim_t = int32_t;
    using bond_idx_t = int32_t;
    using bond_label_t = int32_t;
    using ten_size_t = size_t;
    using elem_t = float;
    using elem_coor_t = int32_t;
    using elem_coors_t = std::vector<int32_t>;
    using real_t = float;
    using real_ten_t = gqten::tensor<float>;
    using cplx_t = std::complex<float>;
    using cplx_ten_t = gqten::tensor<std::complex<float>>;
    using context_handle_t = gqten_handle;
  };

  // traits for gqten::tensor<double>
  template <>
  struct tensor_traits<gqten::tensor<double>> {
    using ten_t = gqten::tensor<double>;
    using order_t = int32_t;
    using shape_t = std::vector<int32_t>;
    using bond_dim_t = int32_t;
    using bond_idx_t = size_t;
    using bond_label_t = int32_t;
    using ten_size_t = size_t;
    using elem_t = double;
    using elem_coor_t = int32_t;
    using elem_coors_t = std::vector<int32_t>;
    using real_t = double;
    using real_ten_t = gqten::tensor<double>;
    using cplx_t = std::complex<double>;
    using cplx_ten_t = gqten::tensor<std::complex<double>>;
    using context_handle_t = gqten_handle;
  };

  // traits for gqten::tensor<std::complex<float>>
  template <>
  struct tensor_traits<gqten::tensor<std::complex<float>>> {
    using ten_t = gqten::tensor<std::complex<float>>;
    using order_t = int32_t;
    using shape_t = std::vector<int32_t>;
    using bond_dim_t = int32_t;
    using bond_idx_t = int32_t;
    using bond_label_t = int32_t;
    using ten_size_t = size_t;
    using elem_t = std::complex<float>;
    using elem_coor_t = int32_t;
    using elem_coors_t = std::vector<int32_t>;
    using real_t = float;
    using real_ten_t = gqten::tensor<float>;
    using cplx_t = std::complex<float>;
    using cplx_ten_t = gqten::tensor<std::complex<float>>;
    using context_handle_t = gqten_handle;
  };

  // traits for gqten::tensor<std::complex<double>>
  template <>
  struct tensor_traits<gqten::tensor<std::complex<double>>> {
    using ten_t = gqten::tensor<std::complex<double>>;
    using order_t = int32_t;
    using shape_t = std::vector<int32_t>;
    using bond_dim_t = int32_t;
    using bond_idx_t = int32_t;
    using bond_label_t = int32_t;
    using ten_size_t = size_t;
    using elem_t = std::complex<double>;
    using elem_coor_t = int32_t;
    using elem_coors_t = std::vector<int32_t>;
    using real_t = double;
    using real_ten_t = gqten::tensor<double>;
    using cplx_t = std::complex<double>;
    using cplx_ten_t = gqten::tensor<std::complex<double>>;
    using context_handle_t = gqten_handle;
  };

} // namespace tci

#endif // TCI_GQTEN_TRAITS_H
