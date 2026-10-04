// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <random>
#include <chrono>

#include "tcapi/tcapi.h"

#include "option.h"

/**
Main function to text singular value decomposition for a random tensor
*/
int main(int argc, char * argv[]) {

  /**
     Setup calculation
   */
  Option option = generate_option(argc,argv);

  ContextHandle ctx;
  tcapi::create_context(ctx);

  ContextHandleReal ctr;
  tcapi::create_context(ctr);

  std::mt19937 engine(option.seed);
  Shape shape = option.shape;
  Tensor A;
  A = tcapi::random<Tensor>(ctx,shape,engine);
  Tensor U;
  Tensor V;
  RealTensor S;
  Real trunc_err;

  auto time_start_svd = std::chrono::high_resolution_clock::now();
  
  tcapi::trunc_svd(ctx,A,option.num_rows,
		 U,S,V,trunc_err,
		 option.chi_max,
		 option.s_min);

  auto time_end_svd = std::chrono::high_resolution_clock::now();
  auto elapsed_count_svd = std::chrono::duration_cast<std::chrono::microseconds>(time_end_svd-time_start_svd).count();
  std::cout << " Execution time for SVD = " << elapsed_count_svd << " (ms) " << std::endl;
  
  std::cout << " Bond Dimension = " << tcapi::shape(ctr,S)[0] << std::endl;

  tcapi::for_each(ctr,S,[](auto & elem){ elem = std::sqrt(elem); });
  Tensor Z;
  tcapi::convert(ctr,S,ctx,Z);
  tcapi::show(ctx,Z);

  return 0;
}
