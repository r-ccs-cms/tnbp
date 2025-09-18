/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/constructor.h
@brief header to define TCI constructors and destructors for gqten::tensor<ElemT>
*/
#ifndef TCI_GQTEN_CONSTRUCTOR_H
#define TCI_GQTEN_CONSTRUCTOR_H

#include <vector>
#include <complex>
#include <utility>     // std::move
#include <functional>  // std::invoke
#include <iterator>    // std::next
#include <type_traits> // std::is_nothrow_*

namespace tci {

//====================== allocate / zeros / eye / fill ====================

template <typename ElemT>
inline void allocate(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape,
  gqten::tensor<ElemT>& a) {
  a = gqten::tensor<ElemT>(shape);
}

template <typename ElemT>
inline gqten::tensor<ElemT> allocate(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape) {
  return gqten::tensor<ElemT>(shape);
}

template <typename ElemT>
inline void zeros(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape,
  gqten::tensor<ElemT>& a) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<ElemT> data(size, ElemT{});
  a = gqten::tensor<ElemT>(shape, data);
}

template <typename ElemT>
inline gqten::tensor<ElemT> zeros(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<ElemT> data(size, ElemT{});
  return gqten::tensor<ElemT>(shape, data);
}

template <typename ElemT>
inline void eye(
  context_handle_t<gqten::tensor<ElemT>>&,
  const bond_dim_t<gqten::tensor<ElemT>> N,
  gqten::tensor<ElemT>& a) {
  const size_t NN = static_cast<size_t>(N);
  std::vector<ElemT> data(NN*NN, ElemT{});
  for (size_t k = 0; k < NN; ++k) data[k + NN*k] = ElemT{1};
  a = gqten::tensor<ElemT>({N, N}, data);
}

template <typename ElemT>
inline gqten::tensor<ElemT> eye(
  context_handle_t<gqten::tensor<ElemT>>&,
  const bond_dim_t<gqten::tensor<ElemT>> N) {
  const size_t NN = static_cast<size_t>(N);
  std::vector<ElemT> data(NN*NN, ElemT{});
  for (size_t k = 0; k < NN; ++k) data[k + NN*k] = ElemT{1};
  return gqten::tensor<ElemT>({N, N}, data);
}

template <typename ElemT>
inline void fill(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape,
  const elem_t<gqten::tensor<ElemT>> v,
  gqten::tensor<ElemT>& a) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<ElemT> data(size, static_cast<ElemT>(v));
  a = gqten::tensor<ElemT>(shape, data);
}

template <typename ElemT>
inline gqten::tensor<ElemT> fill(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape,
  const elem_t<gqten::tensor<ElemT>> v) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<ElemT> data(size, static_cast<ElemT>(v));
  return gqten::tensor<ElemT>(shape, data);
}

//================== assign_from_container (2 形) ==================

/**
 * coors2idx: (const std::vector<size_t>&) -> size_t
 *  address = coor[0] + shape[0]*coor[1] + shape[1]*shape[0]*coor[2] + ...
 */
template <typename ElemT, typename RandomIt, typename Func>
inline void assign_from_container(
  context_handle_t<gqten::tensor<ElemT>>&,
  const shape_t<gqten::tensor<ElemT>>& shape,
  RandomIt init_elems_begin,
  Func&& coors2idx,
  gqten::tensor<ElemT>& a) {
  if (shape.empty()) {
    // For scaler 
    a = gqten::tensor<ElemT>({}, init_elems_begin);
    return;
  }

  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<ElemT>  data(size, ElemT{});
  std::vector<size_t> coor(shape.size());

  for (size_t i = 0; i < size; ++i) {
    size_t tmp = i;
    for (size_t k = 0; k < shape.size(); ++k) {
      coor[k] = tmp % shape[k];
      tmp    /= shape[k];
    }
    const size_t address =
      std::invoke(std::forward<Func>(coors2idx), coor);

    data[i] = *(init_elems_begin + static_cast<std::ptrdiff_t>(address));
    // for arbitrary iterator
    // data[i] = *std::next(init_elems_begin, static_cast<std::ptrdiff_t>(address));
  }

  a = gqten::tensor<ElemT>(shape, data);
}

template <typename ElemT, typename RandomIt, typename Func>
inline gqten::tensor<ElemT> assign_from_container(
  context_handle_t<gqten::tensor<ElemT>>& ctx,
  const shape_t<gqten::tensor<ElemT>>&    shape,
  RandomIt                                 init_elems_begin,
  Func&&                                    coors2idx) {
  gqten::tensor<ElemT> a;
  assign_from_container(ctx, shape, init_elems_begin,
                        std::forward<Func>(coors2idx), a);
  return a;
}

//========================== random (2 形) =========================

template <typename ElemT, typename RandNumGen>
inline void random(
  context_handle_t<gqten::tensor<ElemT>> & ctx,
  const shape_t<gqten::tensor<ElemT>>& shape,
  RandNumGen& gen,
  gqten::tensor<ElemT>& a) {
  a = gqten::tensor<ElemT>(shape);
  a.Random(3, {}, gen);
}

template <typename ElemT, typename RandNumGen>
inline gqten::tensor<ElemT> random(
  context_handle_t<gqten::tensor<ElemT>> & ctx,
  const shape_t<gqten::tensor<ElemT>> & shape,
  RandNumGen & gen) {
  gqten::tensor<ElemT> a;
  random(ctx, shape, gen, a);
  return a;
}

//================== copy / move / clear ==================

template <typename ElemT>
inline void copy(
  context_handle_t<gqten::tensor<ElemT>>&,
  const gqten::tensor<ElemT>& orig,
  gqten::tensor<ElemT>& dist) {
  dist = orig;
}

template <typename ElemT>
inline gqten::tensor<ElemT> copy(
  context_handle_t<gqten::tensor<ElemT>>&,
  const gqten::tensor<ElemT>& orig) {
  return orig;
}

template <typename ElemT>
inline void move(
  context_handle_t<gqten::tensor<ElemT>>&,
  gqten::tensor<ElemT>& from,
  gqten::tensor<ElemT>& to) noexcept(
  std::is_nothrow_move_assignable_v<gqten::tensor<ElemT>>) {
  to = std::move(from);
}

template <typename ElemT>
inline gqten::tensor<ElemT> move(
  context_handle_t<gqten::tensor<ElemT>>&,
  gqten::tensor<ElemT>& from) noexcept(
    std::is_nothrow_move_constructible_v<gqten::tensor<ElemT>>) {
  return std::move(from);
}

template <typename ElemT>
inline void clear(
  context_handle_t<gqten::tensor<ElemT>>&,
  gqten::tensor<ElemT>& a) {
  a = gqten::tensor<ElemT>();
}

} // namespace tci

#endif // TCI_GQTEN_CONSTRUCTOR_H
