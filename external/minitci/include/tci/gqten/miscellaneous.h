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
  template <typename ElemT>
  void create_context(
       context_handle_t<gqten::tensor<ElemT>> &ctx) {}

  
  /**
  template <typename ContextHandleT>
  void destroy_context(ContextHandleT &ctx);
  */
  template <typename ElemT>
  void destroy_context(
       context_handle_t<gqten::tensor<ElemT>> &ctx) {}

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
  template <typename ElemT,
          typename RandomIt,
          typename Func>
  void to_container(
       context_handle_t<gqten::tensor<ElemT>> &ctx,
       const gqten::tensor<ElemT> & a,
       RandomIt first,
       Func &&coors2idx) {
    auto shape = a.Shape();
    if( shape.empty() ) {
      a.GetScale(*first);
      return;
    }

    size_t size = 1;
    for(auto d : shape) size *= d;

    std::vector<size_t> coor(shape.size());

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
  template <typename ElemT>
  void show(
       context_handle_t<gqten::tensor<ElemT>> & ctx,
       const gqten::tensor<ElemT> & a) {
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
  template <typename ElemT>
  bool eq(
       context_handle_t<gqten::tensor<ElemT>> & ctx,
       const gqten::tensor<ElemT> &a,
       const gqten::tensor<ElemT> &b,
       const real_t<gqten::tensor<ElemT>> epsilon) {
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
    using RealT = typename tensor_traits<gqten::tensor<ElemT>>::real_t;
    RealT diff = 0.0;
    const ElemT * raw_a = a.GetRaw();
    const ElemT * raw_b = b.GetRaw();
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
  template <typename Elem1T, typename Elem2T>
  void convert(
       context_handle_t<gqten::tensor<Elem1T>> & ctx1,
       const gqten::tensor<Elem1T> & t1,
       context_handle_t<gqten::tensor<Elem2T>> & ctx2,
       gqten::tensor<Elem2T> & t2) {
    const std::vector<size_t> shape = t1.Shape();
    const Elem1T * raw1 = t1.GetRaw();
    const size_t size = t1.Size();
    Elem2T * raw2 = static_cast<Elem2T*> malloc(sizeof(Elem2T)*size);
    for(size_t i=0; i < size; ++i) {
      raw2[i] = static_cast<Elem2T>(gqten::GetReal(raw1[i]));
    }
    t2 = gqten::tensor<Elem2T>(shape,raw2);
  }
  
  
}
#endif
