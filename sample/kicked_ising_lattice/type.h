#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <numbers>

#ifdef USE_SINGLE
using Tensor = typename gqten::tensor<std::complex<float>>;
#else
using Tensor = typename gqten::tensor<std::complex<double>>;
#endif

using Elem = typename tci::tensor_traits<Tensor>::elem_t;
using Real = typename tci::tensor_traits<Tensor>::real_t;
using ContextHandle = typename tci::tensor_traits<Tensor>::context_handle_t;
using BondDim = typename tci::tensor_traits<Tensor>::bond_dim_t;

template <typename T>
inline constexpr T PI_v = std::numbers::pi_v<T>;
