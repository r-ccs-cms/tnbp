#ifndef TCI_GQTEN_MISCELLANEOUS_H
#define TCI_GQTEN_MISCELLANEOUS_H

// C++17-compatible rewrite of miscellaneous routines for gqten backend.
// - replaces `requires ...` with SFINAE on template parameters
// - preserves original interfaces and semantics
// - depends on tci/traits.h providing: is_gqten_handle_v<T>, is_gqten_tensor_v<T>, 
//   context_handle_t<T>, elem_t<T>, real_t<T>, bond_dim_t<T>.

#include <vector>
#include <functional>  // std::invoke
#include <type_traits>
#include <cstddef>
#include <cstdlib>
#include <cmath>       // std::abs

namespace tci {

  //-------------------------------------------------------------------------
  // create_context / destroy_context
  //-------------------------------------------------------------------------
  template <typename ContextHandleT,
            typename std::enable_if<is_gqten_handle_v<ContextHandleT>, int>::type = 0>
  inline void create_context(ContextHandleT &ctx) {
    ctx.text = std::string("created");
  }

  template <typename ContextHandleT,
            typename std::enable_if<is_gqten_handle_v<ContextHandleT>, int>::type = 0>
  inline void destroy_context(ContextHandleT &ctx) {
    ctx.text = std::string("destroied");
  }

  //-------------------------------------------------------------------------
  // to_container
  //   coors2idx: (const std::vector<bond_dim_t<TenT>>&)->size_t
  //-------------------------------------------------------------------------
  template <typename TenT, typename RandomIt, typename Func,
            typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void to_container(
       context_handle_t<TenT> &/*ctx*/,
       const TenT & a,
       RandomIt first,
       Func &&coors2idx) {
    auto shape = a.Shape();
    if (shape.empty()) {
      *first = a.GetScale();
      return;
    }

    size_t size = 1;
    for (auto d : shape) size *= d;

    std::vector<bond_dim_t<TenT>> coor(shape.size());

    for (size_t i = 0; i < size; ++i) {
      size_t tmp = i;
      for (size_t k = 0; k < shape.size(); ++k) {
        coor[k] = tmp % shape[k];
        tmp    /= shape[k];
      }
      const size_t address =
        std::invoke(std::forward<Func>(coors2idx), coor);

      *(first + static_cast<std::ptrdiff_t>(address)) = a.GetElem(coor);
    }
  }

  //-------------------------------------------------------------------------
  // show
  //-------------------------------------------------------------------------
  template <typename TenT,
            typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline void show(
       context_handle_t<TenT> & /*ctx*/,
       const TenT & a) {
    a.FormattedPrint();
  }

  //-------------------------------------------------------------------------
  // eq
  //-------------------------------------------------------------------------
  template <typename TenT,
            typename std::enable_if<is_gqten_tensor_v<TenT>, int>::type = 0>
  inline bool eq(
       context_handle_t<TenT> & /*ctx*/,
       const TenT &a,
       const TenT &b,
       const real_t<TenT> epsilon) {
    auto shape_a = a.Shape();
    auto shape_b = b.Shape();
    if (shape_a.size() != shape_b.size()) {
      return false;
    }
    auto rank = a.Rank();
    for (size_t k = 0; k < rank; ++k) {
      if (shape_a[k] != shape_b[k]) return false;
    }
    using RealT = typename tensor_traits<TenT>::real_t;
    RealT diff = RealT(0);
    const elem_t<TenT> * raw_a = a.GetRaw();
    const elem_t<TenT> * raw_b = b.GetRaw();
    size_t size  = a.Size();
    for (size_t i = 0; i < size; i++) {
      const auto d = raw_a[i] - raw_b[i];
      diff += std::abs(d * d);
    }
    if (diff > epsilon) return false;
    return true;
  }

  //-------------------------------------------------------------------------
  // convert
  //-------------------------------------------------------------------------
  template <typename Ten1T, typename Ten2T,
            typename std::enable_if<is_gqten_tensor_v<Ten1T> && is_gqten_tensor_v<Ten2T>, int>::type = 0>
  inline void convert(
       context_handle_t<Ten1T> & /*ctx1*/,
       const Ten1T & t1,
       context_handle_t<Ten2T> & /*ctx2*/,
       Ten2T & t2) {
    const std::vector<bond_dim_t<Ten2T>> shape = t1.Shape();
    const elem_t<Ten1T> * raw1 = t1.GetRaw();
    const size_t size = t1.Size();
    elem_t<Ten2T> * raw2 = static_cast<elem_t<Ten2T>*>(std::malloc(sizeof(elem_t<Ten2T>) * size));
    for (size_t i = 0; i < size; ++i) {
      raw2[i] = static_cast<elem_t<Ten2T>>(gqten::GetReal(raw1[i]));
    }
    t2 = Ten2T(shape, raw2);
  }

} // namespace tci

#endif // TCI_GQTEN_MISCELLANEOUS_H
