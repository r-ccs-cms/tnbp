// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "qasm/any.h"
#include "qasm/utility.h"

#include "option.h"

#ifdef USE_COMPLEX
using Tensor = typename gqten::tensor<gqten::Cplx_16>;
#else
using Tensor = typename gqten::tensor<gqten::Real_8>;
#endif

using Elem = typename tci::tensor_traits<Tensor>::elem_t;
using Real = typename tci::tensor_traits<Tensor>::real_t;
using ContextHandle = typename tci::tensor_traits<Tensor>::context_handle_t;

int main(int argc, char * argv[]) {

  int mpi_err = MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_master = 0;
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  
  Option options =generate_option(argc,argv);

  std::string qasm_source;
  if( mpi_rank == mpi_master ) {
    qasm_source = qasm::read_all_file(option.circuit);
  }
  MpiBcast(qasmsource,mpi_master,comm);
  qasm::Program qasm_program = qasm::parse_any(qasm_source);

  ContextHandle ctx;
  tci::create_context(ctx);
  std::vector<std::pair<int,int>> edges;
  std::vector<std::vector<std::pair<int,int>>> layer_edges;
  
  if( option.device == std::string("ibm_kobe") ) {
    edges = tnbp::bond_ibm_kobe();
    layer_edges = tnbp::parallel_bond_ibm_kobe();
  }
  
  auto TPO = tnbp::QasmToTPO(ctx,qasm_program,edges,option.num_gates);

  
  

  return 0;
}
