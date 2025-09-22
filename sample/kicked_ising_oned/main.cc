// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>

#include "qasm/any.h"
#include "qasm/utility.h"
#include "pauli/sparse_pauli.h"
#include "tnbp/tnbp.h"

#include "type.h"
#include "option.h"
#include "circuit.h"
#include "measop.h"
#include "timestamp.h"

int main(int argc, char * argv[]) {

  int mpi_err = MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_master = 0;
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  
  Option options =generate_options(argc,argv);

  if( mpi_rank == 0 ) {
    std::cout << "### Input options " << std::endl;
    cout_options(options);
  }

  ContextHandle ctx;
  tci::create_context<Tensor>(ctx);
  std::vector<std::pair<int,int>> edges;
  std::vector<std::vector<std::pair<int,int>>> layer_edges;
  
  edges = tnbp::bond_oned_lattice(options.num_qubits,0);
  layer_edges = tnbp::parallel_bond_oned_lattice(options.num_qubits,0);

  if( mpi_rank == 0 ) {
    std::cout << " " << make_timestamp() << " Start construction of TPO " << std::endl;
  }
  std::vector<Tensor> TPO = CircuitTPO<Tensor>(ctx,options.num_qubits,
			    options.Jz,options.hz,options.hx,options.dt,comm);

  std::vector<int> qubit = tnbp::GetSiteIndexFromBond(edges);
  std::vector<Tensor> V;
  std::vector<int> SiteIdx;
  std::vector<int> Site_To_MpiRank;
  std::vector<Tensor> E;
  std::vector<int> EdgeIdx;
  std::vector<int> pdim(qubit.size(),2);
  std::vector<Tensor> F;

  if( mpi_rank == 0 ) {
    std::cout << " " << make_timestamp() << " start construction of Initial TPS " << std::endl;
  }
  tnbp::InitTensorProductState(ctx,edges,pdim,
			       V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm);
  
  Real tolerance;
  for(uint32_t step=0; step < options.num_steps; step++) {
    if( mpi_rank == 0 ) {
      std::cout << " " << make_timestamp() << " start AbsorbTPO " << std::endl;
    }
    tnbp::AbsorbTPO(
	  ctx,TPO,edges,
	  V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,comm);

    if( mpi_rank == 0 ) {
      std::cout << " " << make_timestamp() << " start belief propagation loop " << std::endl;
    }
    for(uint32_t k=0; k < options.max_bp_iterations; k++) {
      
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp() << " start " << k << "th belief propagation " << std::endl;
      }
      for(size_t p=0; p < layer_edges.size(); p++) {
	tnbp::BeliefPropagation(
	    ctx,edges,
	    V,SiteIdx,Site_To_MpiRank,
	    layer_edges[p],E,EdgeIdx,comm,F);
      }
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp() << " update edge tensors " << std::endl;
      }
      auto itE = E.begin();
      for(const auto & et : F) {
	tci::copy(ctx,et,*itE++);
      }
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp() << " calculate condition " << std::endl;
      }
      tnbp::BeliefPropagationCondition(
	    ctx,edges,V,SiteIdx,
	    Site_To_MpiRank,E,EdgeIdx,
	    comm,tolerance);
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp() << " error from belief propagation condition = "
		  << tolerance << std::endl;
      }
      if( tolerance < options.bp_tolerance ) {
	break;
      }
    }
    tnbp::Truncation(ctx,edges,
		     V,SiteIdx,Site_To_MpiRank,
		     E,EdgeIdx,comm,
		     options.max_bond_dim,
		     options.sv_min,
		     options.truncation_error);
    
  }

  MPI_Finalize();
  
  return 0;
}
