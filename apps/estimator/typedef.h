#pragma once
#include <complex>
#include "tcapi/tcapi.h"

#if defined(USE_TCAPI_CUDA) && defined(USE_CYTNX)
#error "Choose only one tensor backend: USE_TCAPI_CUDA or USE_CYTNX"
#endif

#ifdef USE_TCAPI_CUDA
#ifdef USE_SINGLE
using Tensor = tcapi::cuda::Tensor<std::complex<float>>;
#else
using Tensor = tcapi::cuda::Tensor<std::complex<double>>;
#endif
#elif defined(USE_CYTNX)
using Tensor = tcapi::CytnxTensor<cytnx::cytnx_complex128>;
#else
#ifdef USE_SINGLE
using Tensor = typename gqten::tensor<std::complex<float>>;
#else
using Tensor = typename gqten::tensor<std::complex<double>>;
#endif
#endif

using Elem = typename tcapi::tensor_traits<Tensor>::elem_t;
using Real = typename tcapi::tensor_traits<Tensor>::real_t;
using ContextHandle = typename tcapi::tensor_traits<Tensor>::context_handle_t;
using BondDim = typename tcapi::tensor_traits<Tensor>::bond_dim_t;

inline float GetReal(float a) { return a; }
inline double GetReal(double a) { return a; }
inline float GetReal(std::complex<float> a) { return a.real(); }
inline double GetReal(std::complex<double> a) { return a.real(); }

