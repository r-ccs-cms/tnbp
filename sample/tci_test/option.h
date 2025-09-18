#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef USE_COMPLEX
#ifdef USE_SINGLE
using Tensor = typename gqten::tensor<std::complex<float>>;
#else
using Tensor = typename gqten::tensor<std::complex<double>>;
#endif
#else
#ifdef USE_SINGLE
using Tensor = typename gqten::tensor<float>;
#else
using Tensor = typename gqten::tensor<double>;
#endif
#endif

using ContextHandle = typename tci::tensor_traits<Tensor>::context_handle_t;
using Elem = typename tci::tensor_traits<Tensor>::elem_t;
using Real = typename tci::tensor_traits<Tensor>::real_t;
using Rank = typename tci::tensor_traits<Tensor>::rank_t;
using Shape = typename tci::tensor_traits<Tensor>::shape_t;
using BondDim = typename tci::tensor_traits<Tensor>::bond_dim_t;
using BondLabel = typename tci::tensor_traits<Tensor>::bond_label_t;
using RealTensor = typename tci::tensor_traits<Tensor>::real_ten_t;

struct Option {
  Rank num_rows = 1;
  std::vector<BondDim> shape = {10,2,10};
  BondDim chi_min = 1;
  BondDim chi_max = 10;
  Real trunc_err = 1.0e-8;
  Real s_min = 1.0e-8;
  uint32_t seed = 1729;
};
  
Option generate_option(int argc, char * argv[]) {
  Option option;
  for(int i=0; i < argc; i++) {
    if( std::string(argv[i]) == "--num_rows" ) {
      option.num_rows = static_cast<Rank>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--shape" ) {
      std::stringstream ss(argv[++i]);
      std::string token;
      while ( std::getline(ss, token, ',')) {
	option.shape.push_back(static_cast<BondDim>(std::stoi(token)));
      }
    }
    if( std::string(argv[i]) == "--chi_min" ) {
      option.chi_min = static_cast<BondDim>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--chi_max" ) {
      option.chi_max = static_cast<BondDim>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--trunc_err" ) {
      option.trunc_err = static_cast<Real>(std::atof(argv[++i]));
    }
    if( std::string(argv[i]) == "--s_min" ) {
      option.s_min   = static_cast<Real>(std::atof(argv[++i]));
    }
    if( std::string(argv[i]) == "--seed" ) {
      option.seed = std::atoi(argv[++i]);
    }
  }
  return option;
}
