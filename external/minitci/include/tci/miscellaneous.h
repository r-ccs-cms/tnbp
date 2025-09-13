/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/miscellaneous.h
@brief header to declare miscellaneous routines
 */
#ifndef TCI_MISCELLANEOUS_H
#define TCI_MISCELLANEOUS_H

namespace tci {

  template <typename ContextHandleT>
  void create_context(ContextHandleT &ctx);

  template <typename ContextHandleT>
  void destroy_context(ContextHandleT &ctx);

  template <typename TenT,
          typename RandomIt,
          typename Func>
  void to_container(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       RandomIt first,
       Func &&coors2idx);

  template <typename TenT>
  void show(
       context_handle_t<TenT> &ctx,
       const TenT &a);

  template <typename TenT>
  bool eq(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const TenT &b,
       const elem_t<TenT> epsilon);

  template <typename Ten1T, typename Ten2T>
  void convert(
       context_handle_t<Ten1T> &ctx1,
       const Ten1T &t1,
       context_handle_t<Ten2T> &ctx2,
       Ten2T &t2);
  
}

#endif
