/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/constructor.h
@brief header to declare constructors
 */
#ifndef TCI_CONSTRUCTOR_H
#define TCI_CONSTRUCTOR_H
namespace tci {

  template <typename TenT>
  void allocate(
       context_handle_t<TenT> & ctx,
       const shape_t<TenT> & shape,
       TenT & a);

  template <typename TenT>
  TenT allocate(
      context_handle_t<TenT> & ctx,
      const shape_t<TenT> & shape);

  template <typename TenT>
  void zeros(
       context_handle_t<TenT> & ctx,
       const shape_t<TenT> & shape,
       TenT & a);

  template <typename TenT>
  TenT zeros(
       context_handle_t<TenT> & ctx,
       const shape_t<TenT> & shape);

  template <typename TenT,
	    typename RandomIt,
	    typename Func>
  void assign_from_container(
       context_handle_t<TenT> &ctx,
       const shape_t<TenT> &shape,
       RandomIt init_elems_begin,
       Func &&coors2idx,
       TenT &a);

  template <typename TenT,
	    typename RandomIt,
	    typename Func>
  TenT assign_from_container(
        context_handle_t<TenT> &ctx,
	const shape_t<TenT> &shape,
	RandomIt init_elems_begin,
	Func &&coors2idx);

  template <typename TenT,
	    typename RandNumGen>
  void random(
       context_handle_t<TenT> &ctx,
       const shape_t<TenT> &shape,
       RandNumGen &gen,
       TenT &a);

  template <typename TenT,
	    typename RandNumGen>
  TenT random(
       context_handle_t<TenT> &ctx,
       const shape_t<TenT> &shape,
       RandNumGen &gen);

  template <typename TenT>
  void eye(
       context_handle_t<TenT> &ctx,
       const bond_dim_t<TenT> N,
       TenT &a);

  template <typename TenT>
  TenT eye(
       context_handle_t<TenT> &ctx,
       const bond_dim_t<TenT> N);

  template <typename TenT>
  void fill(
       context_handle_t<TenT> &ctx,
       const shape_t<TenT> &shape,
       const elem_t<TenT> v,
       TenT &a);

  template <typename TenT>
  TenT fill(
       context_handle_t<TenT> &ctx,
       const shape_t<TenT> &shape,
       const elem_t<TenT> v);

  template <typename TenT>
  void copy(
       context_handle_t<TenT> &ctx,
       const TenT &orig,
       TenT &dist);

  template <typename TenT>
  TenT copy(
       context_handle_t<TenT> &ctx,
       const TenT &orig);

  template <typename TenT>
  void move(
       context_handle_t<TenT> &ctx,
       TenT &from,
       TenT &to);

  template <typename TenT>
  TenT move(
       context_handle_t<TenT> &ctx,
       TenT &from);

  template <typename TenT>
  void clear(
       context_handle_t<TenT> &ctx,
       TenT &a);
  
}
#endif
