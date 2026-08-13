#ifdef USE_CYTNX
using Tensor = tci::CytnxTensor<cytnx::cytnx_complex128>;
#else
using Tensor = typename gqten::tensor<std::complex<double>>;
#endif

using Elem = typename tci::tensor_traits<Tensor>::elem_t;
using Real = typename tci::tensor_traits<Tensor>::real_t;
using BondDim = typename tci::tensor_traits<Tensor>::bond_dim_t;
using ContextHandle = typename tci::tensor_traits<Tensor>::context_handle_t;
