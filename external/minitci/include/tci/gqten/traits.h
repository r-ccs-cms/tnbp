/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/traits.h
@brief header to define traits for gqten::tensor<ElemT>
*/
#ifndef TCI_GQTEN_TRAITS_H
#define TCI_GQTEN_TRAITS_H

namespace tci {

  template<class T> struct is_gqten_tensor : std::false_type {};
  template<class E> struct is_gqten_tensor<gqten::tensor<E>> : std::true_type {};
  template<class T>
  inline constexpr bool is_gqten_tensor_v = is_gqten_tensor<T>::value;
  
  struct gqten_handle {
    std::string text = std::string("none");
  };
  
  // traits for gqten::tensor<float>
  template <>
  struct tensor_traits<gqten::tensor<float>> {
    using ten_t = gqten::tensor<float>;
    using rank_t = int32_t;
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
    using rank_t = int32_t;
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
    using rank_t = int32_t;
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
    using rank_t = int32_t;
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
  
}

#endif
