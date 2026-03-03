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
  tci::create_context(ctx);
  std::vector<std::pair<int,int>> edges;
  std::vector<std::vector<std::pair<int,int>>> layer_edges;

  if( options.lattice == std::string("oned") ) {
    edges = tnbp::bond_oned_lattice(options.L,0);
    layer_edges = tnbp::parallel_bond_oned_lattice(options.L,0);
  } else if ( options.lattice == std::string("honeycomb") ) {
    edges = tnbp::bond_honeycomb_lattice(options.L,options.L);
    layer_edges = tnbp::parallel_bond_honeycomb_lattice(options.L,options.L);
  } else if ( options.lattice == std::string("square") ) {
    edges = tnbp::bond_square_lattice(options.L,options.L);
    layer_edges = tnbp::GreedyEdgeLayering(edges);
  }

  if( mpi_rank == 0 ) {
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
    std::cout << " " << make_timestamp() << " Start construction of TPO " << std::endl;
    
  }
  std::vector<Tensor> TPO = CircuitTPO<Tensor>(ctx,edges,
			    options.Jz,options.hz,options.hx,options.dt,comm);

  std::vector<int> qubit = tnbp::GetSiteIndexFromBond(edges);
  std::vector<Tensor> V;
  std::vector<int> SiteIdx;
  std::map<int,int> Site_To_MpiRank;
  std::vector<Tensor> E;
  std::vector<int> EdgeIdx;
  std::vector<int> pdim(qubit.size(),2);
  std::vector<Tensor> F;
  std::vector<BondDim> res_bond_dim;
  std::vector<Real> res_truncation_error;

  if( mpi_rank == 0 ) {
    std::cout << " " << make_timestamp() << " start construction of Initial TPS " << std::endl;
  }
  tnbp::InitTensorProductState(ctx,edges,pdim,
			       V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm);

  if ( mpi_rank == 0 ) {
    std::cout << " " << make_timestamp()
	      << " start construction of meas operators " << std::endl;
  }

  std::vector<int> meas_site;
  std::vector<Tensor> Ox;
  std::vector<Tensor> Oz;
  std::vector<std::pair<int,int>> meas_edge;
  std::vector<Tensor> Oj;

  measops(ctx,edges,meas_site,Ox,Oz,meas_edge,Oj);

  std::cout.precision(16);
  Real tolerance;
  for(uint32_t step=0; step < options.num_steps; step++) {
    if( mpi_rank == 0 ) {
      std::cout << " " << make_timestamp()
		<< " start AbsorbTPO " << std::endl;
    }
    tnbp::AbsorbTPO(
	  ctx,TPO,edges,
	  V,SiteIdx,Site_To_MpiRank,E,EdgeIdx,comm);

    if( mpi_rank == 0 ) {
      std::cout << " " << make_timestamp()
		<< " start belief propagation loop " << std::endl;
    }
    for(uint32_t k=0; k < options.max_bp_iterations; k++) {
      
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp()
		  << " start " << k
		  << "th belief propagation " << std::endl;
      }
      for(size_t p=0; p < layer_edges.size(); p++) {
	tnbp::BeliefPropagation(
	    ctx,edges,
	    V,SiteIdx,Site_To_MpiRank,
	    layer_edges[p],E,EdgeIdx,comm,F);
      }
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp()
		  << " update edge tensors " << std::endl;
      }
      auto itE = E.begin();
      for(const auto & et : F) {
	tci::copy(ctx,et,*itE++);
      }
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp()
		  << " calculate condition " << std::endl;
      }
      tnbp::BeliefPropagationCondition(
	    ctx,edges,V,SiteIdx,
	    Site_To_MpiRank,E,EdgeIdx,
	    comm,tolerance);
      if( mpi_rank == 0 ) {
	std::cout << " " << make_timestamp()
		  << " error from belief propagation condition at time step "
		  << step << " = "
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
		     options.truncation_error,
		     res_bond_dim,
		     res_truncation_error);
    for(size_t k=0; k < EdgeIdx.size(); k++) {
      std::cout << " " << make_timestamp()
		<< " truncation error for edge ("
		<< edges[EdgeIdx[k]].first << ","
		<< edges[EdgeIdx[k]].second
		<< ") = " << res_truncation_error[k]
		<< " bond dim = "
		<< res_bond_dim[k] << std::endl;
    }

    std::cout << " " << make_timestamp()
	      << " start measurements " << std::endl;
    
    auto res_x = tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm,
			       meas_site,Ox);
    
    for(size_t i=0; i < meas_site.size(); i++) {
      if( Site_To_MpiRank[meas_site[i]] == mpi_rank ) {
	std::cout << " " << make_timestamp()
		  << " measurement of X at site " << meas_site[i]
		  << " time step " << step
		  << " = " << GetReal(res_x[i]) << std::endl;
      }
    }
    
    auto res_z = tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm,
			       meas_site,Oz);
    
    for(size_t i=0; i < meas_site.size(); i++) {
      if( Site_To_MpiRank[meas_site[i]] == mpi_rank ) {
	std::cout << " " << make_timestamp()
		  << " measurement of Z at site " << meas_site[i]
		  << " time step " << step
		  << " = " << GetReal(res_z[i]) << std::endl;
      }
    }
    
    auto res_j = tnbp::Measure(ctx,edges,V,SiteIdx,Site_To_MpiRank,
			       E,EdgeIdx,comm,
			       meas_edge,Oj);
    
    for(size_t m=0; m < meas_edge.size(); m++) {
      int site_a = edges[m].first;
      int site_b = edges[m].second;
      int mpi_rank_a = Site_To_MpiRank[site_a];
      int mpi_rank_b = Site_To_MpiRank[site_b];
      int mpi_type = 0;
      if ( mpi_rank_a == mpi_rank ) {
	mpi_type += 2;
      }
      if ( mpi_rank_b == mpi_rank ) {
	mpi_type += 1;
      }
      if ( mpi_type > 1 ) {
	std::cout << " " << make_timestamp()
		  << " measurement of ZZ at sites (" << meas_edge[m].first
		  << "," << meas_edge[m].second
		  << ") time step " << step
		  << " = " << GetReal(res_j[m]) << std::endl;
      }
    }
    
  }


  MPI_Finalize();
  
  return 0;
}
