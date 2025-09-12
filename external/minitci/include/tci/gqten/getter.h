/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/getter.h
@brief header to define TCI getters for gqten::tensor<ElemT>
 */
#ifndef TCI_GQTEN_GETTER_H
#define TCI_GQTEN_GETTER_H

namespace tci {

  /**
     template <typename TenT>
     rank_t<TenT> rank(
         context_handle_t<TenT> & ctx,
	 const TenT & a);
   */
  template <>
  rank_t<gqten::tensor<float>> rank(
	 context_handle_t<gqten::tensor<float>> & ctx,
	 const gqten::tensor<float> & a) {
    return a.Rank();
  }

  template <>
  rank_t<gqten::tensor<double>> rank(
	 context_handle_t<gqten::tensor<double>> & ctx,
	 const gqten::tensor<double> & a) {
    return a.Rank();
  }
  
  template <>
  rank_t<gqten::tensor<std::complex<float>>> rank(
	 context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
	 const gqten::tensor<std::complex<float>> & a) {
    return a.Rank();
  }

  template <>
  rank_t<gqten::tensor<std::complex<double>>> rank(
	 context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
	 const gqten::tensor<std::complex<double>> & a) {
    return a.Rank();
  }

  /**
     template <typename TenT>
     shape_t<TenT> shape(
     context_handle_t<TenT> & ctx,
     const TenT & a);
   */

  template <>
  shape_t<gqten::tensor<float>> shape(
       context_handle_t<gqten::tensor<float>> & ctx,
       const gqten::tensor<float> & a) {
    return a.Shape();
  }
		      
  template <>
  shape_t<gqten::tensor<double>> shape(
       context_handle_t<gqten::tensor<double>> & ctx,
       const gqten::tensor<double> & a) {
    return a.Shape();
  }
		      
  template <>
  shape_t<gqten::tensor<std::complex<float>>> shape(
       context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
       const gqten::tensor<std::complex<float>> & a) {
    return a.Shape();
  }
		      
  template <>
  shape_t<gqten::tensor<std::complex<double>>> shape(
       context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
       const gqten::tensor<std::complex<double>> & a) {
    return a.Shape();
  }

  /**
     template <typename TenT>
     ten_size_t<TenT> size(
     context_handle_t<TenT> & ctx,
     const TenT & a);
   */
  template <>
  ten_size_t<gqten::tensor<float>> size(
       context_handle_t<gqten::tensor<float>> & ctx,
       const gqten::tensor<float> & a) {
    return a.Size();
  }
  
  template <>
  ten_size_t<gqten::tensor<double>> size(
       context_handle_t<gqten::tensor<double>> & ctx,
       const gqten::tensor<double> & a) {
    return a.Size();
  }
  
  template <>
  ten_size_t<gqten::tensor<std::complex<float>>> size(
       context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
       const gqten::tensor<std::complex<float>> & a) {
    return a.Size();
  }
  
  template <>
  ten_size_t<gqten::tensor<std::complex<double>>> size(
       context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
       const gqten::tensor<std::complex<double>> & a) {
    return a.Size();
  }
  
  /**
     template <typename TenT>
     std::size_t size_bytes(
     context_handle_t<TenT> & ctx,
     const TenT & a);
   */
  template <>
  std::size_t size_bytes(
	context_handle_t<gqten::tensor<float>> & ctx,
	const gqten::tensor<float> & a) {
    return a.Size()*sizeof(float);
  }
  
  template <>
  std::size_t size_bytes(
	context_handle_t<gqten::tensor<double>> & ctx,
	const gqten::tensor<double> & a) {
    return a.Size()*sizeof(double);
  }
  
  template <>
  std::size_t size_bytes(
	context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
	const gqten::tensor<std::complex<float>> & a) {
    return a.Size()*sizeof(std::complex<float>);
  }
  
  template <>
  std::size_t size_bytes(
	context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
	const gqten::tensor<std::complex<double>> & a) {
    return a.Size()*sizeof(std::complex<double>);
  }

  /**
     template <typename TenT>
     void get_elem(
     context_handle_t<TenT> & ctx,
     const TenT & a,
     const elem_coors_t<TenT> & coors,
     elem_t<TenT> & elem);
   */
  template <>
  void get_elem(
	context_handle_t<gqten::tensor<float>> & ctx,
	const gqten::tensor<float> & a,
	const elem_coors_t<gqten::tensor<float>> & coors,
	elem_t<gqten::tensor<float>> & elem) {
    elem = a.GetElem(coors);
  }
  
  template <>
  void get_elem(
	context_handle_t<gqten::tensor<double>> & ctx,
	const gqten::tensor<double> & a,
	const elem_coors_t<gqten::tensor<double>> & coors,
	elem_t<gqten::tensor<double>> & elem) {
    elem = a.GetElem(coors);
  }
  
  template <>
  void get_elem(
	context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
	const gqten::tensor<std::complex<float>> & a,
	const elem_coors_t<gqten::tensor<std::complex<float>>> & coors,
	elem_t<gqten::tensor<std::complex<float>>> & elem) {
    elem = a.GetElem(coors);
  }
  
  template <>
  void get_elem(
	context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
	const gqten::tensor<std::complex<double>> & a,
	const elem_coors_t<gqten::tensor<std::complex<double>>> & coors,
	elem_t<gqten::tensor<std::complex<double>>> & elem) {
    elem = a.GetElem(coors);
  }

  /**
     template <typename TenT>
     elem_t<TenT> get_elem(
     context_handle_t<TenT> & ctx,
     const TenT & a,
     const elem_coors_t<TenT> & coors);
   */
  template <>
  elem_t<gqten::tensor<float>> get_elem(
	context_handle_t<gqten::tensor<float>> & ctx,
	const gqten::tensor<float> & a,
	const elem_coors_t<gqten::tensor<float>> & coors) {
    return a.GetElem(coors);
  }

  template <>
  elem_t<gqten::tensor<double>> get_elem(
	context_handle_t<gqten::tensor<double>> & ctx,
	const gqten::tensor<double> & a,
	const elem_coors_t<gqten::tensor<double>> & coors) {
    return a.GetElem(coors);
  }
  
  template <>
  elem_t<gqten::tensor<std::complex<float>>> get_elem(
	context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
	const gqten::tensor<std::complex<float>> & a,
	const elem_coors_t<gqten::tensor<std::complex<float>>> & coors) {
    return a.GetElem(coors);
  }

  template <>
  elem_t<gqten::tensor<std::complex<double>>> get_elem(
	context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
	const gqten::tensor<std::complex<double>> & a,
	const elem_coors_t<gqten::tensor<std::complex<double>>> & coors) {
    return a.GetElem(coors);
  }


}

#endif
