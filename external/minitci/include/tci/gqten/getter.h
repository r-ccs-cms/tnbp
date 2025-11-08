#ifndef TCI_GQTEN_GETTER_H
#define TCI_GQTEN_GETTER_H

// C++17-compatible rewrite of getters for gqten backend.
// - replaces `requires is_gqten_tensor_v<TenT>` with SFINAE on template parameter
// - preserves original interfaces and semantics
// - depends on tci/traits.h providing order_t<T>, shape_t<T>, ten_size_t<T>, elem_t<T>, elem_coors_t<T>, context_handle_t<T>,
//   and is_gqten_tensor_v<T>.

#include <cstddef>   // std::size_t
#include <complex>   // std::complex
#include <type_traits>

namespace tci {

  // order
  template <typename TenT>
  inline order_t<TenT> order(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return a.Rank();
  }

  // shape
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline shape_t<TenT> shape(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return a.Shape();
  }

  // size (number of elements)
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline ten_size_t<TenT> size(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
    return a.Size();
  }

  // size in bytes
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline std::size_t size_bytes(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a) {
        return static_cast<std::size_t>(a.Size()) * sizeof(elem_t<TenT>);
  }

  // get_elem (out-param)
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void get_elem(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a,
      const elem_coors_t<TenT>& coors,
      elem_t<TenT>& elem) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    ShapeT shape_a = a.Shape();
    if (shape_a.size() == 0) {
      elem = a.GetScale();
    } else {
      elem = a.GetElem(coors);
    }
  }

  // get_elem (by value)
  template <typename TenT, typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline elem_t<TenT> get_elem(
      context_handle_t<TenT>& /*ctx*/,
      const TenT& a,
      const elem_coors_t<TenT>& coors) {
    using ElemT = typename tci::tensor_traits<TenT>::elem_t;
    using ShapeT = typename tci::tensor_traits<TenT>::shape_t;
    ShapeT shape_a = a.Shape();
    if (shape_a.size() == 0) {
      ElemT scale = a.GetScale();
      return scale;
    }
    return a.GetElem(coors);
  }

} // namespace tci

#endif // TCI_GQTEN_GETTER_H
