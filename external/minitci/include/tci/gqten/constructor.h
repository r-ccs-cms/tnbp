/// This file is a part of r-ccs-cms/tnbp
/**
@file tci/gqten/constructor.h
@brief header to define TCI constructors and destructors for gqten::tensor<ElemT>
 */
#ifndef TCI_GQTEN_CONSTRUCTOR_H
#define TCI_GQTEN_CONSTRUCTOR_H

namespace tci {

  /**
     template <typename TenT>
     void allocate(
     context_handle_t<TenT> & ctx,
     const shape_t<TenT> & shape,
     TenT & a);
   */
  template <>
  void allocate(
       context_handle_t<gqten::tensor<float>> & ctx,
       const shape_t<gqten::tensor<float>> & shape,
       gqten::tensor<float> & a) {
    a = gqten::tensor<float>(shape);
  }
  
  template <>
  void allocate(
       context_handle_t<gqten::tensor<double>> & ctx,
       const shape_t<gqten::tensor<double>> & shape,
       gqten::tensor<double> & a) {
    a = gqten::tensor<double>(shape);
  }
  
  template <>
  void allocate(
       context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
       const shape_t<gqten::tensor<std::complex<float>>> & shape,
       gqten::tensor<std::complex<float>> & a) {
    a = gqten::tensor<std::complex<float>>(shape);
  }
  
  template <>
  void allocate(
       context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
       const shape_t<gqten::tensor<std::complex<double>>> & shape,
       gqten::tensor<std::complex<double>> & a) {
    a = gqten::tensor<std::complex<double>>(shape);
  }

  /**
     template <typename TenT>
     TenT allocate(
     context_handle_t<TenT> & ctx,
     const shape_t<TenT> & shape);
   */
  template <>
  gqten::tensor<float> allocate(
	context_handle_t<gqten::tensor<float>> & ctx,
	const shape_t<gqten::tensor<float>> & shape) {
    return gqten::tensor<float>(shape);
  }
  
  template <>
  gqten::tensor<double> allocate(
	context_handle_t<gqten::tensor<double>> & ctx,
	const shape_t<gqten::tensor<double>> & shape) {
    return gqten::tensor<double>(shape);
  }
  
  template <>
  gqten::tensor<std::complex<float>> allocate(
	context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
	const shape_t<gqten::tensor<std::complex<float>>> & shape) {
    return gqten::tensor<std::complex<float>>(shape);
  }
  
  template <>
  gqten::tensor<std::complex<double>> allocate(
	context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
	const shape_t<gqten::tensor<std::complex<double>>> & shape) {
    return gqten::tensor<std::complex<double>>(shape);
  }

  /**
     template <typename TenT>
     void zeros(
     context_handle_t<TenT> & ctx,
     const shape_t<TenT> & shape,
     TenT & a);
   */
  template <>
  void zeros(
       context_handle_t<gqten::tensor<float>> & ctx,
       const shape_t<gqten::tensor<float>> & shape,
       gqten::tensor<float> & a) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<float> data(size,0.0);
    a = gqten::tensor<float>(shape,data.data());
  }
  
  template <>
  void zeros(
       context_handle_t<gqten::tensor<double>> & ctx,
       const shape_t<gqten::tensor<double>> & shape,
       gqten::tensor<double> & a) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<double> data(size,0.0);
    a = gqten::tensor<double>(shape,data.data());
  }
  
  template <>
  void zeros(
       context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
       const shape_t<gqten::tensor<std::complex<float>>> & shape,
       gqten::tensor<std::complex<float>> & a) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<std::complex<float>> data(size,static_cast<std::complex<float>>(0.0));
    a = gqten::tensor<std::complex<float>>(shape,data.data());
  }
  
  template <>
  void zeros(
       context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
       const shape_t<gqten::tensor<std::complex<double>>> & shape,
       gqten::tensor<std::complex<double>> & a) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<std::complex<double>> data(size,static_cast<std::complex<double>>(0.0));
    a = gqten::tensor<std::complex<double>>(shape,data.data());
  }

  /**
     template <typename TenT>
     TenT zeros(
     context_handle_t<TenT> & ctx,
     const shape_t<TenT> & shape);
   */
  template <>
  gqten::tensor<float> zeros(
       context_handle_t<gqten::tensor<float>> & ctx,
       const shape_t<gqten::tensor<float>> & shape) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<float> data(size,static_cast<float>(0.0));
    return gqten::tensor<float>(shape,data.data());
  }

  template <>
  gqten::tensor<double> zeros(
       context_handle_t<gqten::tensor<double>> & ctx,
       const shape_t<gqten::tensor<double>> & shape) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<double> data(size,static_cast<double>(0.0));
    return gqten::tensor<double>(shape,data.data());
  }
  
  template <>
  gqten::tensor<std::complex<float>> zeros(
       context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
       const shape_t<gqten::tensor<std::complex<float>>> & shape) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<std::complex<float>> data(size,static_cast<std::complex<float>>(0.0));
    return gqten::tensor<std::complex<float>>(shape,data.data());
  }

  template <>
  gqten::tensor<std::complex<double>> zeros(
       context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
       const shape_t<gqten::tensor<std::complex<double>>> & shape) {
    size_t size = 1;
    for(const auto & dim : shape) {
      size *= dim;
    }
    std::vector<std::complex<double>> data(size,static_cast<std::complex<double>>(0.0));
    return gqten::tensor<std::complex<double>>(shape,data.data());
  }

  

  
}

#endif
