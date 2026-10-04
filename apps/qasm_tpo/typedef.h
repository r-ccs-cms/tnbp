#ifdef USE_CYTNX
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

