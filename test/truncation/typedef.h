#ifdef USE_CYTNX
using Tensor = tcapi::CytnxTensor<cytnx::cytnx_complex128>;
#else
using Tensor = typename gqten::tensor<std::complex<double>>;
#endif

using Elem = typename tcapi::tensor_traits<Tensor>::elem_t;
using Real = typename tcapi::tensor_traits<Tensor>::real_t;
using BondDim = typename tcapi::tensor_traits<Tensor>::bond_dim_t;
using ContextHandle = typename tcapi::tensor_traits<Tensor>::context_handle_t;
