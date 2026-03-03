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
  
  Option options =generate_options(argc,argv);

  if( mpi_rank == mpi_master ) {
    cout_options(options);
  }

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
  
  if( options.backend == std::string("ibm_kobe") ) {
    edges = tnbp::bond_ibm_kobe();
    layer_edges = tnbp::parallel_bond_ibm_kobe();
  } else if ( options.backend == std::string("default") ) {
    edges = tnbp::EdgesFromQasm(qasm_program);
    layer_edges = tnbp::GreedyEdgeLayering(edges);
  }

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
  if( options.num_gates.size() != 0 ) {
    TPO = tnbp::QasmToTPO<Tensor>(ctx,qasm_program,edges,
				  options.num_gates);
    if( mpi_rank == mpi_master ) {
      std::cout << " " << make_timestamp()
		<< " manual gate-number layering of QasmToTPO was performed " << std::endl;
    }
  } else {
    TPO = tnbp::QasmToTPO<Tensor>(ctx,qasm_program,edges);
    if( mpi_rank == mpi_master ) {
      std::cout << " " << make_timestamp()
		<< " barrier-based layering of QasmToTPO was performed " << std::endl;
    }
  }

  if( options.do_opt_tpo ) {
    for(size_t l=0; l < TPO.size(); l++) {
      tnbp::OptTPObySVD(ctx,edges,TPO[l],options.eps_opt_tpo);
      std::cout << " " << make_timestamp()
		<< " optimization of tpo for layer " << l << " is performed " << std::endl;
    }
  }

  std::vector<int> qubit = tnbp::GetSiteIndexFromBond(edges);
  std::vector<Tensor> V;
  std::vector<int> SiteIdx;
  std::map<int,int> Site_To_MpiRank;
  std::vector<Tensor> E;
  std::vector<int> EdgeIdx;
  std::vector<int> pdim(qubit.size(),2);
  std::vector<Tensor> F;

  tnbp::InitTensorProductState(ctx,edges,pdim,
			       V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm);
  

  Real tolerance;
  std::vector<BondDim> res_bond_dim;
  std::vector<Real> res_truncation_error;

  if( mpi_rank == mpi_master ) {
    std::cout << " " << make_timestamp()
	      << " Start circuit contractions " << std::endl;
  }
  
  for(size_t t=0; t < TPO.size(); t++) {
    /**
       Attach operator tensor to site tensor
     */
    tnbp::AbsorbTPO(
	  ctx,TPO[t],edges,
	  V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,
	  comm);

    /**
       Perform belief propagation loop
     */
    for(size_t k=0; k < options.max_bp_iterations; k++) {
      
      for(size_t p=0; p < layer_edges.size(); p++) {
	tnbp::BeliefPropagation(
	      ctx,edges,V,SiteIdx,
	      Site_To_MpiRank,
	      layer_edges[p],E,EdgeIdx,
	      comm,F);
      }
      auto itE = E.begin();
      for(const auto & et : F) {
	tci::copy(ctx,et,*itE++);
      }
      tnbp::BeliefPropagationCondition(
	    ctx,edges,V,SiteIdx,
	    Site_To_MpiRank,E,EdgeIdx,
	    comm,tolerance);
      if( mpi_rank == mpi_master ) {
	std::cout << " " << make_timestamp()
		  << " error from belief propagation condition at step "
		  << t << " = "
		  << tolerance << std::endl;
      }
      if( tolerance < options.bp_tolerance ) {
	break;
      }
    }
    
    /**
       Perform truncation
    */
    tnbp::Truncation(ctx,edges,
		     V,SiteIdx,Site_To_MpiRank,
		     E,EdgeIdx,comm,
		     options.max_bond_dim,
		     options.sv_min,
		     options.truncation_error,
		     res_bond_dim,
		       res_truncation_error);
    for(size_t k=0; k < EdgeIdx.size(); k++) {
      std::cout << " " << make_timestamp()
		<< " truncation error for edge ("
		<< edges[EdgeIdx[k]].first << ","
		<< edges[EdgeIdx[k]].second << ") = "
		<< res_truncation_error[k] << " bond dim = "
		<< res_bond_dim[k] << std::endl;
    }
  }

  /**
     Measurement
   */

  if( mpi_rank == mpi_master ) {
    std::cout << " " << make_timestamp()
	      << " start measurements " << std::endl;
  }
  auto spo = pauli::load_sparse_pauli_op<Elem>(options.sparsepauli);
  std::vector<Tensor> exOp;
  std::vector<std::vector<int>> exSite;
  tnbp::SparsePauliToTensorOp(ctx,spo,exOp,exSite);
  if( mpi_rank == mpi_master ) {
    std::cout << " " << make_timestamp()
	      << " number of measrement operators = " << exOp.size() << std::endl;
  }
  std::vector<Tensor> exOp_one;
  std::vector<int> exSite_one;
  std::vector<Tensor> exOp_two;
  std::vector<std::pair<int,int>> exPair_two;
  tnbp::ClassifyOps(ctx,edges,exOp,exSite,
	      exOp_one,exSite_one,
	      exOp_two,exPair_two);

  if( mpi_rank == mpi_master ) {
    std::cout << " " << make_timestamp()
	      << " generated single site operator " << std::endl;
  }

  std::vector<Elem> exVal_one =
    tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
		  E,EdgeIdx,comm,exSite_one,exOp_one);

  std::cout.precision(16);
  for(int site_address=0; site_address < exSite_one.size(); site_address++) {
    if( mpi_rank == Site_To_MpiRank[exSite_one[site_address]] ) {
      std::cout << " Expectation value of single site operator at site "
		<< exSite_one[site_address]
		<< " " << GetReal(exVal_one[site_address]) << std::endl;
    }
  }

  std::vector<Elem> exVal_two_rank =
    tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
		  E,EdgeIdx,comm,exPair_two,exOp_two);
  std::vector<Elem> exVal_two(exPair_two.size());
  
  MPI_Datatype DataT = tnbp::GetMpiType<Elem>();
  MPI_Allreduce(exVal_two_rank.data(),exVal_two.data(),
		exVal_two.size(),DataT,MPI_SUM,comm);

  if( mpi_rank == mpi_master ) {
    for(int edge_address=0; edge_address < exPair_two.size(); edge_address++) {
      std::cout << " Expectation value of two site operator at sites ("
		<< exPair_two[edge_address].first
		<< "," << exPair_two[edge_address].second
		<< ") "
		<< GetReal(exVal_two[edge_address]) << std::endl;
    }
  }

  MPI_Finalize();
  
  return 0;
}
