/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/miscellaneous.h
@brief header to define miscellaneous routines for gqten::tensor<ElemT>
 */
#ifndef TCI_GQTEN_MISCELLANEOUS_H
#define TCI_GQTEN_MISCELLANEOUS_H

#include <cstdlib>

namespace tci {
  /**
  template <typename ContextHandleT>
  void create_context(ContextHandleT &ctx);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  inline void create_context(
    context_handle_t<TenT> &ctx) {
    ctx.text = std::string("created");
  }

  
  /**
  template <typename ContextHandleT>
  void destroy_context(ContextHandleT &ctx);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void destroy_context(
    context_handle_t<TenT> &ctx) {
    ctx.text = std::string("destroied");
  }

  /**
  template <typename TenT,
          typename RandomIt,
          typename Func>
  void to_container(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       RandomIt first,
       Func &&coors2idx);
  */
  template <typename TenT,
          typename RandomIt,
          typename Func>
  requires is_gqten_tensor_v<TenT>
  void to_container(
       context_handle_t<TenT> &ctx,
       const TenT & a,
       RandomIt first,
       Func &&coors2idx) {
    auto shape = a.Shape();
    if( shape.empty() ) {
      *first = a.GetScale();
      return;
    }

    size_t size = 1;
    for(auto d : shape) size *= d;

    std::vector<bond_dim_t<TenT>> coor(shape.size());

    for(size_t i=0; i < size; ++i) {
      size_t tmp = i;
      for(size_t k=0; k < shape.size(); ++k) {
	coor[k] = tmp % shape[k];
	tmp    /= shape[k];
      }
      const size_t address =
	std::invoke(std::forward<Func>(coors2idx),coor);
      
      *(first + static_cast<std::ptrdiff_t>(address)) =
	a.GetElem(coor);
    }
  }
  
  /**
  template <typename TenT>
  void show(
       context_handle_t<TenT> &ctx,
       const TenT &a);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  void show(
       context_handle_t<TenT> & ctx,
       const TenT & a) {
    a.FormattedPrint();
  }

  /**
  template <typename TenT>
  bool eq(
       context_handle_t<TenT> &ctx,
       const TenT &a,
       const TenT &b,
       const real_t<TenT> epsilon);
  */
  template <typename TenT>
  requires is_gqten_tensor_v<TenT>
  bool eq(
       context_handle_t<TenT> & ctx,
       const TenT &a,
       const TenT &b,
       const real_t<TenT> epsilon) {
    auto shape_a = a.Shape();
    auto shape_b = b.Shape();
    if( shape_a.size() != shape_b.size() ) {
      return false;
    }
    bool res = true;
    auto rank = a.Rank();
    for(size_t k=0; k < rank; ++k) {
      if( shape_a[k] != shape_b[k] ) {
	return false;
      }
    }
    using RealT = typename tensor_traits<TenT>::real_t;
    RealT diff = 0.0;
    const elem_t<TenT> * raw_a = a.GetRaw();
    const elem_t<TenT> * raw_b = b.GetRaw();
    size_t size  = a.Size();
    for(size_t i=0; i < size; i++) {
      diff += std::abs((raw_a[i]-raw_b[i])*(raw_a[i]-raw_b[i]));
    }
    if( diff > epsilon ) {
      return false;
    }
    return true;
  }
       
  /**
  template <typename Ten1T, typename Ten2T>
  void convert(
       context_handle_t<Ten1T> &ctx1,
       const Ten1T &t1,
       context_handle_t<Ten2T> &ctx2,
       Ten2T &t2);
  */
  template <typename Ten1T, typename Ten2T>
  requires is_gqten_tensor_v<Ten1T> && is_gqten_tensor_v<Ten2T>
  void convert(
       context_handle_t<Ten1T> & ctx1,
       const Ten1T & t1,
       context_handle_t<Ten2T> & ctx2,
       Ten2T & t2) {
    const std::vector<bond_dim_t<Ten2T>> shape = t1.Shape();
    const elem_t<Ten1T> * raw1 = t1.GetRaw();
    const size_t size = t1.Size();
    elem_t<Ten2T> * raw2 = static_cast<elem_t<Ten2T>*>(malloc(sizeof(elem_t<Ten2T>)*size));
    for(size_t i=0; i < size; ++i) {
      raw2[i] = static_cast<elem_t<Ten2T>>(gqten::GetReal(raw1[i]));
    }
    t2 = Ten2T(shape,raw2);
  }
  
  
}
#endif
