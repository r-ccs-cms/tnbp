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
