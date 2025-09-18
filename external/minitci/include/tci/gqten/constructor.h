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

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void allocate(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape,
  TenT& a) {
  a = TenT(shape);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline TenT allocate(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape) {
  return TenT(shape);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void zeros(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape,
  TenT& a) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<elem_t<TenT>> data(size, elem_t<TenT>{});
  a = TenT(shape, data);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline TenT zeros(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<elem_t<TenT>> data(size, elem_t<TenT>{});
  return TenT(shape, data);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void eye(
  context_handle_t<TenT>&,
  const bond_dim_t<TenT> N,
  TenT& a) {
  const size_t NN = static_cast<size_t>(N);
  std::vector<elem_t<TenT>> data(NN*NN, elem_t<TenT>{});
  for (size_t k = 0; k < NN; ++k) data[k + NN*k] = elem_t<TenT>{1};
  a = TenT({N, N}, data);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline TenT eye(
  context_handle_t<TenT>&,
  const bond_dim_t<TenT> N) {
  const size_t NN = static_cast<size_t>(N);
  std::vector<elem_t<TenT>> data(NN*NN, elem_t<TenT>{});
  for (size_t k = 0; k < NN; ++k) data[k + NN*k] = elem_t<TenT>{1};
  return TenT({N, N}, data);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void fill(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape,
  const elem_t<TenT> v,
  TenT& a) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<elem_t<TenT>> data(size, static_cast<elem_t<TenT>>(v));
  a = TenT(shape, data);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline TenT fill(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape,
  const elem_t<TenT> v) {
  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<elem_t<TenT>> data(size, static_cast<elem_t<TenT>>(v));
  return TenT(shape, data);
}

//================== assign_from_container (2 形) ==================

/**
 * coors2idx: (const std::vector<size_t>&) -> size_t
 *  address = coor[0] + shape[0]*coor[1] + shape[1]*shape[0]*coor[2] + ...
 */
template <typename TenT, typename RandomIt, typename Func>
requires is_gqten_tensor_v<TenT>
inline void assign_from_container(
  context_handle_t<TenT>&,
  const shape_t<TenT>& shape,
  RandomIt init_elems_begin,
  Func&& coors2idx,
  TenT& a) {
  if (shape.empty()) {
    // For scaler
    elem_t<TenT> * data = static_cast<elem_t<TenT>*>(
		    std::malloc(static_cast<size_t>(1)));
    *data = *init_elems_begin;
    a = TenT(shape, data);
    return;
  }

  size_t size = 1;
  for (auto d : shape) size *= d;

  std::vector<elem_t<TenT>>  data(size, elem_t<TenT>{});
  std::vector<bond_dim_t<TenT>> coor(shape.size());

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

  a = TenT(shape, data);
}

template <typename TenT, typename RandomIt, typename Func>
requires is_gqten_tensor_v<TenT>
inline TenT assign_from_container(
  context_handle_t<TenT>& ctx,
  const shape_t<TenT>&    shape,
  RandomIt                                 init_elems_begin,
  Func&&                                    coors2idx) {
  TenT a;
  assign_from_container(ctx, shape, init_elems_begin,
                        std::forward<Func>(coors2idx), a);
  return a;
}

//========================== random (2 形) =========================

template <typename TenT, typename RandNumGen>
requires is_gqten_tensor_v<TenT>
inline void random(
  context_handle_t<TenT> & ctx,
  const shape_t<TenT>& shape,
  RandNumGen& gen,
  TenT & a) {
  a = TenT(shape);
  a.Random(3, {}, gen);
}

template <typename TenT, typename RandNumGen>
requires is_gqten_tensor_v<TenT>
inline TenT random(
  context_handle_t<TenT> & ctx,
  const shape_t<TenT> & shape,
  RandNumGen & gen) {
  TenT a;
  random(ctx, shape, gen, a);
  return a;
}

//================== copy / move / clear ==================

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void copy(
  context_handle_t<TenT>&,
  const TenT& orig,
  TenT& dist) {
  dist = orig;
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline TenT copy(
  context_handle_t<TenT>&,
  const TenT& orig) {
  return orig;
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void move(
  context_handle_t<TenT>&,
  TenT& from,
  TenT& to) noexcept(
  std::is_nothrow_move_assignable_v<TenT>) {
  to = std::move(from);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline TenT move(
  context_handle_t<TenT>&,
  TenT& from) noexcept(
    std::is_nothrow_move_constructible_v<TenT>) {
  return std::move(from);
}

template <typename TenT>
requires is_gqten_tensor_v<TenT>
inline void clear(
  context_handle_t<TenT>&,
  TenT& a) {
  a = TenT();
}

} // namespace tci

#endif // TCI_GQTEN_CONSTRUCTOR_H
