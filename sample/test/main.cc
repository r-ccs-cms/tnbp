#include <random>
#include <chrono>

#include "tci.h"
#include "tci_gqten.h"

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
using Shape = typename tci::tensor_traits<Tensor>::shape_t;


int main(int argc, char * argv[]) {

  ContextHandle ctx;
  tci::create_context(ctx);

  uint32_t seed = 1729;
  std::mt19937 engine(seed);
  Shape shape = {3,4,6};
  Tensor A = tci::random<Tensor>(ctx,shape,engine);
  
  return 0;
}
