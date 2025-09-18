#include <complex>
#include <vector>
#include <string>
#include "gqten/gqten.h"


namespace tci {

  struct gqten_handle {
    std::string text;
  };

  template <>
  struct tensor_traits<gqten::tensor<float>> {
    using context_handle_t = gqten_handle;
    using shape_t = std::vector<int>;
  };

  template <>
  struct tensor_traits<gqten::tensor<double>> {
    using context_handle_t = gqten_handle;
    using shape_t = std::vector<int>;
  };
  
  template <>
  struct tensor_traits<gqten::tensor<std::complex<float>>> {
    using context_handle_t = gqten_handle;
    using shape_t = std::vector<int>;
  };
  
  template <>
  struct tensor_traits<gqten::tensor<std::complex<double>>> {
    using context_handle_t = gqten_handle;
    using shape_t = std::vector<int>;
  };


  inline void create_context(gqten_handle & ctx) {
    ctx.text = std::string("created");
  }

  template <typename ElemT, typename RandNumGen>
  void random(
      gqten_handle & ctx,
      const shape_t<gqten::tensor<ElemT>> & shape,
      RandNumGen & gen,
      gqten::tensor<ElemT> & a) {
    a = gqten::tensor<ElemT>(shape);
    a.Random(3, {}, gen);
  }

  template <typename TenT, typename RandNumGen,
	    
	    >
  gqten::tensor<float> random(
      context_handle_t<gqten::tensor<float>> & ctx,
      const shape_t<gqten::tensor<float>> & shape,
      RandNumGen & gen) {
    gqten::tensor<float> a = gqten::tensor<float>(shape);
    a.Random(3, {}, gen);
    return a;
  }

  template <typename RandNumGen>
  gqten::tensor<double> random(
      context_handle_t<gqten::tensor<double>> & ctx,
      const shape_t<gqten::tensor<double>> & shape,
      RandNumGen & gen) {
    gqten::tensor<double> a = gqten::tensor<double>(shape);
    a.Random(3, {}, gen);
    return a;
  }

  template <typename RandNumGen>
  gqten::tensor<std::complex<float>> random(
      context_handle_t<gqten::tensor<std::complex<float>>> & ctx,
      const shape_t<gqten::tensor<std::complex<float>>> & shape,
      RandNumGen & gen) {
    gqten::tensor<std::complex<float>> a = gqten::tensor<std::complex<float>>(shape);
    a.Random(3, {}, gen);
    return a;
  }

  template <typename RandNumGen>
  gqten::tensor<std::complex<double>> random(
      context_handle_t<gqten::tensor<std::complex<double>>> & ctx,
      const shape_t<gqten::tensor<std::complex<double>>> & shape,
      RandNumGen & gen) {
    gqten::tensor<std::complex<double>> a = gqten::tensor<std::complex<double>>(shape);
    a.Random(3, {}, gen);
    return a;
  }

  
}
