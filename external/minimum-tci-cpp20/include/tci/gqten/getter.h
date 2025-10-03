/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/getter.h
@brief header to define TCI getters for TenT
*/
#ifndef TCI_GQTEN_GETTER_H
#define TCI_GQTEN_GETTER_H

#include <cstddef>   // std::size_t
#include <complex>   // std::complex

namespace tci {

  // rank
  template <typename TenT>
  inline rank_t<TenT> rank(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return a.Rank();
  }

  // shape
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  inline shape_t<TenT> shape(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return a.Shape();
  }

  // size (number of elements)
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  inline ten_size_t<TenT> size(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return a.Size();
  }

  // size in bytes
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  inline std::size_t size_bytes(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return static_cast<std::size_t>(a.Size()) * sizeof(elem_t<TenT>);
  }

  // get_elem (out-param)
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  inline void get_elem(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a,
      const elem_coors_t<TenT>& coors,
      elem_t<TenT>& elem) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    ShapeT shape_a = a.Shape();
    if( shape_a.size() == 0 ) {
      elem = a.GetScale();
    } else {
      elem = a.GetElem(coors);
    }
  }

  // get_elem (by value)
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  inline elem_t<TenT> get_elem(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a,
      const elem_coors_t<TenT>& coors) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    ShapeT shape_a = a.Shape();
    if( shape_a.size() == 0 ) {
      ElemT scale = a.GetScale();
      return scale;
    }
    return a.GetElem(coors);
  }

} // namespace tci

#endif // TCI_GQTEN_GETTER_H
