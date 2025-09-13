/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/getter.h
@brief header to define TCI getters for gqten::tensor<ElemT>
*/
#ifndef TCI_GQTEN_GETTER_H
#define TCI_GQTEN_GETTER_H

#include <cstddef>   // std::size_t
#include <complex>   // std::complex

namespace tci {

  // rank
  template <typename ElemT>
  inline rank_t<gqten::tensor<ElemT>> rank(
      context_handle_t<gqten::tensor<ElemT>>& /*ctx*/,
      const gqten::tensor<ElemT>& a) {
    return a.Rank();
  }

  // shape
  template <typename ElemT>
  inline shape_t<gqten::tensor<ElemT>> shape(
      context_handle_t<gqten::tensor<ElemT>>& /*ctx*/,
      const gqten::tensor<ElemT>& a) {
    return a.Shape();
  }

  // size (number of elements)
  template <typename ElemT>
  inline ten_size_t<gqten::tensor<ElemT>> size(
      context_handle_t<gqten::tensor<ElemT>>& /*ctx*/,
      const gqten::tensor<ElemT>& a) {
    return a.Size();
  }

  // size in bytes
  template <typename ElemT>
  inline std::size_t size_bytes(
      context_handle_t<gqten::tensor<ElemT>>& /*ctx*/,
      const gqten::tensor<ElemT>& a) {
    return static_cast<std::size_t>(a.Size()) * sizeof(ElemT);
  }

  // get_elem (out-param)
  template <typename ElemT>
  inline void get_elem(
      context_handle_t<gqten::tensor<ElemT>>& /*ctx*/,
      const gqten::tensor<ElemT>& a,
      const elem_coors_t<gqten::tensor<ElemT>>& coors,
      elem_t<gqten::tensor<ElemT>>& elem) {
    elem = a.GetElem(coors);
  }

  // get_elem (by value)
  template <typename ElemT>
  inline elem_t<gqten::tensor<ElemT>> get_elem(
      context_handle_t<gqten::tensor<ElemT>>& /*ctx*/,
      const gqten::tensor<ElemT>& a,
      const elem_coors_t<gqten::tensor<ElemT>>& coors) {
    return a.GetElem(coors);
  }

} // namespace tci

#endif // TCI_GQTEN_GETTER_H
