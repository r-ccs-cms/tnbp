// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "qasm/any.h"
#include "qasm/utility.h"
#include "pauli/sparse_pauli.h"
#include "tnbp/tnbp.h"

#include "typedef.h"
#include "timestamp.h"
#include "option.h"

int main(int argc, char * argv[]) {

  int mpi_err = MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_master = 0;
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  
  Option options = generate_options(argc,argv);

  std::string qasm_source;
  if( mpi_rank == mpi_master ) {
    qasm_source = qasm::read_all_file(options.circuit);
  }
  tnbp::MpiBcast(qasm_source,mpi_master,comm);
  qasm::Program qasm_program = qasm::parse_any(qasm_source);

  ContextHandle ctx;
  tci::create_context(ctx);
  std::vector<std::pair<int,int>> edges;
  std::vector<std::vector<std::pair<int,int>>> layer_edges;
  
  edges = tnbp::EdgesFromQasm(qasm_program);
  layer_edges = tnbp::GreedyEdgeLayering(edges);

  if( mpi_rank == mpi_master ) {
    std::cout << " " << make_timestamp()
	      << " edges:";
    for(const auto & [u,v] : edges) {
      std::cout << " (" << u << "," << v << ")";
    }
    std::cout << std::endl;
    int m=0;
    for(const auto & layer : layer_edges) {
      std::cout << " " << make_timestamp()
		<< " layer " << m << " edges:";
      for(const auto & [u,v] : layer) {
	std::cout << " (" << u << "," << v << ")";
      }
      std::cout << std::endl;
      m++;
    }
    
    std::cout << " " << make_timestamp()
	      << " Start TPO construction for quantum circuit " << std::endl;
  }

  std::vector<std::vector<Tensor>> TPO;

  TPO = tnbp::QasmToTPO<Tensor>(ctx,qasm_program,edges);
  if( mpi_rank == mpi_master ) {
    std::cout << " " << make_timestamp()
	      << " barrier-based layering of QasmToTPO was performed " << std::endl;
  }

  if( options.do_opt_tpo ) {
    for(size_t l=0; l < TPO.size(); l++) {
      tnbp::OptTPObySVD(ctx,edges,TPO[l],options.opt_tpo_eps);
      std::cout << " " << make_timestamp()
		<< " svd-based optimization for tpo is performed for layer " << l << std::endl;
    }
  }

  if( mpi_rank == mpi_master ) {
    for(size_t l=0; l < TPO.size(); l++) {
      for(size_t i=0; i < TPO[l].size(); i++) {
	std::cout << " Layer " << l << " site " << i << " tensor operator:";
	auto shape = tci::shape(ctx,TPO[l][i]);
	for(size_t k=0; k < shape.size(); k++) {
	  std::cout << ((k==0) ? "{" : ",") << shape[k];
	}
	std::cout << "}" << std::endl;
      }
    }
  }

  MPI_Finalize();
  
  return 0;
}
