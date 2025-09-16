// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "qasm/any.h"
#include "qasm/utility.h"
#include "pauli/sparse_pauli.h"

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

  std::vector<int> qubit = tnbp::GetSiteIndexFromBond(edges);
  std::vector<Tensor> V;
  std::vector<int> SiteIdx;
  std::vector<int> Site_To_MpiRank;
  std::vector<Tensor> E;
  std::vector<int> EdgeIdx;
  std::vector<int> pdim(qubit.size(),2);
  std::vector<Tensor> F;

  tnbp::InitTensorProductState(ctx,edges,pdim,
			       V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm);
  

  for(size_t t=0; t < TPO.size(); t++) {
    /**
       Attach operator tensor to site tensor
     */
    tnbp::AbsrbTPO(
	  ctx,TPO[t],edges,
	  V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,
	  comm);

    /**
       Perform belief propagation loop
     */
    for(size_t k=0; k < option.max_bp_iterations; k++) {
      tnbp::BeliefPropagation(
	    ctx,edges,V,SiteIdx,
	    Site_To_MpiRank,
	    layer_edges,E,EdgeIdx,
	    comm,F);
      auto itE = E.begin();
      for(const auto & et : F) {
	tci::copy(ctx,et,*itE++);
      }
      tnbp::BeliefPropagationCondition(
	    ctx,edges,V,SiteIdx,
	    Site_To_MpiRank,E,EdgeIdx,
	    comm,tolerance);
      if( tolerance < option.tolerance ) {
	break;
      }
    }
  }

  /**
     Measurement
   */
  auto ops = pauli::load_sparse_pauli_op<ElemT>(option.sparsepauli);
  
  
  return 0;
}
