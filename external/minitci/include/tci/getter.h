/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/traits.h
@brief header to define traits
 */
#ifndef TCI_TRAITS_H
#define TCI_TRAITS_H
namespace tci {

  template <typename TenT>
  rank_t<TenT> rank(
       context_handle_t<TenT> & ctx,
       const TenT & a);

  template <typename TenT>
  shape_t<TenT> shape(
       context_handle_t<TenT> & ctx,
       const TenT & a);

  template <typename TenT>
  ten_size_t<TenT> size(
       context_handle_t<TenT> & ctx,
       const TenT & a);

  template <typename TenT>
  std::size_t size_bytes(
       context_handle_t<TenT> & ctx,
       const TenT & a);

  template <typename TenT>
  void get_elem(
       context_handle_t<TenT> & ctx,
       const TenT & a,
       const elem_coors_t<TenT> & coors,
       elem_t<TenT> & elem);

  template <typename TenT>
  elem_t<TenT> get_elem(
       context_handle_t<TenT> & ctx,
       const TenT & a,
       const elem_coors_t<TenT> & coors);

}
#endif
