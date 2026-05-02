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
#include "tnbp/ptns/function/bmps.h"

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
      std::cerr << " input file is not specified " << std::endl;
      MPI_Abort(comm,1);
    }
    edges = read_edges_from_file<int>(options.inputname,true);
    std::cout << " readin edges: size = " << edges.size() << std::endl;
  }
  bcast_edges(edges,0,comm);

  std::vector<Tensor> V;
  std::vector<int> site_idx;
  std::map<int,int> site_to_mpi_rank;

  init_random_tensor(ctx,edges,options.state_bond_dim,
		     V,site_idx,site_to_mpi_rank,
		     options.seed,comm);

  std::vector<std::vector<int>> lines;
  tnbp::bmps_setup_lines(edges,options.init_line,lines);

  std::size_t min_line_idx = 0;
  std::size_t max_line_idx = (lines.size()>0) ? (lines.size()-1) : 0;
  if( mpi_rank == 0 ) {
    for(std::size_t line_idx=min_line_idx; line_idx <= max_line_idx; line_idx++) {
      std::cout << " site in line " << line_idx << ":";
      for(std::size_t site_adrs=0; site_adrs < lines[line_idx].size(); site_adrs++) {
	std::cout << " " << lines[line_idx][site_adrs];
      }
      std::cout << std::endl;
    }
  }

  std::vector<std::pair<int,int>> edges_bmps(edges);
  std::vector<Tensor> T_bmps;
  std::vector<int> site_idx_bmps;
  std::map<int,int> site_to_mpi_rank_bmps(site_to_mpi_rank);
  for(std::size_t line_idx=min_line_idx; line_idx <= max_line_idx; line_idx++) {
    for(std::size_t site_adrs=0; site_adrs < lines[line_idx].size(); site_adrs++) {
      auto site_i = lines[line_idx][site_adrs];
      auto mpi_rank_i = site_to_mpi_rank.at(site_i);
      if( mpi_rank_i == mpi_rank ) {
	auto it_adrs_i = std::find(site_idx.begin(),
				   site_idx.end(),
				   site_i);
	auto adrs_i = std::distance(site_idx.begin(),
				    it_adrs_i);
	T_bmps.push_back(tci::copy(ctx,V[adrs_i]));
	site_idx_bmps.push_back(site_i);
      }
    }
    if( line_idx > min_line_idx ) {
      std::vector<Tensor> E_bmps;
      std::vector<int> edge_idx_bmps;
      std::vector<BondDim> res_bond_dim;
      std::vector<Real> res_trunc_err;
      tnbp::bmps_increment_line(ctx,lines[line_idx-1],lines[line_idx],
				edges_bmps,T_bmps,site_idx_bmps,
				site_to_mpi_rank_bmps,
				comm);
      tnbp::bmps_bp_truncation(ctx,edges_bmps,lines[line_idx],
			       T_bmps,site_idx_bmps,
			       site_to_mpi_rank_bmps,
			       E_bmps,edge_idx_bmps,
			       options.bmps_max_bp_steps,
			       options.bmps_bp_tolerance,
			       options.bmps_bond_dim,
			       options.bmps_sv_min,
			       options.bmps_tg_err,
			       res_bond_dim,
			       res_trunc_err,
			       comm,true);
      // output edges
      if( mpi_rank == 0 ) {
	std::cout << " bmps incrementation from line "
		  << line_idx-1 << " to line " << line_idx
		  << ": edges =";
	for(std::size_t edge_adrs=0;
	    edge_adrs < edges_bmps.size();
	    edge_adrs++) {
	  std::cout << " (" << edges_bmps[edge_adrs].first
		    << "," << edges_bmps[edge_adrs].second
		    << ")";
	}
	std::cout << std::endl;
	if( !options.outputname.empty() ) {
	  auto outputfile = make_step_filename(options.outputname,line_idx);
	  write_edges_to_file(outputfile,edges_bmps,line_idx);
	}
      }
    }
  }
  

  MPI_Finalize();
  return 0;

}
