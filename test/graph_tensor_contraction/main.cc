// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <stdexcept>

#include "mpi.h"

#include "tnbp/tnbp.h"
#include "tnbp/framework/graph.h"
#include "tnbp/framework/helper.h"

#include "typedef.h"
#include "timestamp.h"
#include "option.h"
#include "graphio.h"
#include "gcutils.h"

int main(int argc, char * argv[]) {

  int mpi_err = MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_master = 0;
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);
  
  Option options = generate_options(argc,argv);
  ContextHandle ctx;
  tci::create_context(ctx);
  std::vector<std::pair<int,int>> edges;

  if( mpi_rank == 0 ) {
    cout_options(options);
    if( options.inputname.empty() ) {
      std::cerr << " input filen is not specified " << std::endl;
      MPI_Abort(comm,1);
    }
    edges = read_edges_from_file<int>(options.inputname,true);
    std::cout << " readin edges: size = " << edges.size() << std::endl;
  }
  bcast_edges(edges,0,comm);

  if( edges.size() <= options.num_fin_edges ) {
    std::cerr << " number of final edges > input edges ";
    MPI_Abort(comm,1);
  }

  auto sites = tnbp::GetSiteIndexFromBond(edges);

  // generate tensors compatible to edges
  std::vector<Tensor> W;
  std::vector<int> site_idx;
  std::map<int,int> site_to_mpi_rank;

  init_random_tensor(ctx,edges,options.bond_dim,
		     W,site_idx,site_to_mpi_rank,
		     options.seed_for_ten+mpi_rank,comm);
  
  std::mt19937 rng(options.seed);

  size_t step = 0;
  if( mpi_rank == 0 ) {
    std::cout << " graph at " << step << " step:";
    for(const auto & [u,v] : edges) {
      std::cout << " (" << u << "," << v << ")";
    }
    std::cout << std::endl;
    if( !options.outputname.empty() ) {
      auto outputfile = make_step_filename(options.outputname,step);
      write_edges_to_file(outputfile,edges,step);
    }
  }
  
  while ( edges.size() > options.num_fin_edges ) {
    std::vector<int> cnt_sites(3);
    if( options.schedule_type == std::string("random") ) {
      std::uniform_int_distribution<std::size_t> dist(0,edges.size()-1);
      const auto & edge = edges[dist(rng)];
      cnt_sites[0] = edge.first;
      cnt_sites[1] = edge.second;
      cnt_sites[2] = cnt_sites[0];
      MPI_Bcast(cnt_sites.data(),3,MPI_INT,0,comm);
    }

    if( mpi_rank == 0 ) {
      std::cout << " step " << step
		<< ": contraction between sites (" << cnt_sites[0]
		<< "," << cnt_sites[1] << ")" << std::endl;
    }
    tnbp::graph_tensor_contraction(ctx,edges,W,site_idx,site_to_mpi_rank,
				   cnt_sites[0],cnt_sites[1],cnt_sites[2],comm);
    step++;
    if( mpi_rank == 0 ) {
      if( !options.outputname.empty() ) {
	auto outputfile = make_step_filename(options.outputname,step);
	write_edges_to_file(outputfile,edges,step);
      }
    }
    
    if( mpi_rank == 0 ) {
      std::cout << " graph at " << step << " step:";
      for(const auto & [u,v] : edges) {
	std::cout << " (" << u << "," << v << ")";
      }
      std::cout << std::endl;
    }
    for(int rank=0; rank < mpi_size; rank++) {
      if( mpi_rank == rank ) {
	if( mpi_rank == 0 ) {
	  std::cout << " tensors at " << step << " step:";
	}
	for(size_t i=0; i < site_idx.size(); i++) {
	  if( i != 0 || rank != 0 ) {
	    std::cout << ", ";
	  }
	  auto shape_i = tci::shape(ctx,W[i]);
	  auto order_i = tci::order(ctx,W[i]);
	  for(size_t m=0; m < order_i; m++) {
	    std::cout << ((m==0)? "[" : ",") << shape_i[m];
	  }
	  std::cout << "](site " << site_idx[i] << ")";
	}
      }
      if( rank == mpi_size-1 ) {
	std::cout << std::endl;
      }
      MPI_Barrier(comm);
    }
    
  }

  MPI_Finalize();
  
  return 0;
}
