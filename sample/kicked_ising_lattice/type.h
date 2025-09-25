#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <numbers>

#ifdef USE_CYTNX
using Tensor = typename cytnx::Tensor;
#else
#ifdef USE_SINGLE
using Tensor = typename gqten::tensor<std::complex<float>>;
#else
using Tensor = typename gqten::tensor<std::complex<double>>;
#endif
#endif

using Elem = typename tci::tensor_traits<Tensor>::elem_t;
using Real = typename tci::tensor_traits<Tensor>::real_t;
using ContextHandle = typename tci::tensor_traits<Tensor>::context_handle_t;
using BondDim = typename tci::tensor_traits<Tensor>::bond_dim_t;

template <typename T>
inline constexpr T PI_v = std::numbers::pi_v<T>;

inline float GetReal(float a) { return a; }
inline double GetReal(double a) { return a; }
inline float GetReal(std::complex<float> a) { return a.real(); }
inline double GetReal(std::complex<double> a) { return a.real(); }

