// main.cc
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <stdexcept>
#include <limits>
#include <set>

#include "mpi.h"

#include "tnbp/tnbp.h"
#include "tnbp/framework/graphio.h"
#include "tnbp/ptns/function/mcmc_group_io.h"

#include "typedef.h"
#include "../cuda_device.h"
#include "tnbp/framework/timestamp.h"
#include "option.h"
#include "gcutils.h"
#include "../site_ownership.h"

int main(int argc, char * argv[]) {

  int mpi_err = MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int mpi_master = 0;
  int mpi_size; MPI_Comm_size(comm,&mpi_size);
  int mpi_rank; MPI_Comm_rank(comm,&mpi_rank);

  try {
  Option options = generate_options(argc,argv);
  if(options.help) {
    if(mpi_rank==0) cout_help();
    MPI_Finalize();
    return 0;
  }
  if(mpi_rank==0) cout_options(options);
  std::cout << " " << tnbp::make_timestamp()
            << " start program " << std::endl;
  initialize_tensor_device(comm);
  ContextHandle ctx;
  tcapi::create_context(ctx);
  std::vector<std::pair<int,int>> edges;

  tnbp::PortableCheckpoint portable_descriptor;
  if(!options.portable_input.empty()) {
    if(!options.inputname.empty()) throw std::invalid_argument("choose portable_input or inputname");
    tnbp::portable_detail::phase(comm,[&]{portable_descriptor=tnbp::ReadTPSCheckpoint(options.portable_input);});
    edges=portable_descriptor.edges;
    if(!options.inputedge.empty()) {
      tnbp::portable_detail::phase(comm,[&]{
        if(tnbp::read_edges_from_file<int>(options.inputedge,false)!=edges)
          throw std::invalid_argument("external graph differs from checkpoint graph/order");
      });
    }
  } else {
  if( mpi_rank == 0 ) {
    if( options.inputedge.empty() ) {
      std::cerr << " input file is not specified " << std::endl;
      MPI_Abort(comm,1);
    } else {
      std::cout << " " << tnbp::make_timestamp()
		<< " start readin edges " << std::endl;
    }
    // Checkpoint virtual axes follow file edge order; never sort these edges.
    edges = tnbp::read_edges_from_file<int>(options.inputedge,false);
    std::set<std::pair<int,int>> unique_edges;
    for(const auto& edge:edges)
      if(!unique_edges.insert(edge).second)
        throw std::invalid_argument("duplicate edge in input graph");
    std::cout << " " << tnbp::make_timestamp()
	      << " readin edges: size = " << edges.size()
	      << ", edges =";
    for(const auto & [u,v] : edges) {
      std::cout << " (" << u << "," << v << ")";
    }
    std::cout << "\n";
  }
  tnbp::bcast_edges(edges,0,comm);
  }

  std::size_t num_batches = options.num_batches;

  std::vector<Tensor> W;
  std::vector<int> site_idx;
  std::map<int,int> site_to_mpi_rank;
  std::mt19937 rng(options.seed);
  if(!options.portable_input.empty()) {
    // Independent chains on different owners must not share identical draws.
    // Retain the legacy input/random-state seed behavior for compatibility.
    std::seed_seq sequence{options.seed,static_cast<std::uint32_t>(mpi_rank)};
    rng.seed(sequence);
  }

  if(!options.portable_input.empty()) {
    auto owners=app_ownership::select(portable_descriptor.sites,options.ownerfile,comm,true);
    auto loaded=tnbp::LoadTPSPortable<Tensor>(ctx,options.portable_input,owners,comm,
        tnbp::PortableLoadMode::VerticesOnly);
    W=std::move(loaded.V);site_idx=std::move(loaded.SiteIdx);site_to_mpi_rank=std::move(loaded.owners);
    if(mpi_rank==0) {
      std::cout<<"PORTABLE_INPUT digest="<<loaded.checkpoint.state_digest;
      std::cout<<std::endl;
    }
  } else if( options.inputname.empty() ) {
    site_to_mpi_rank=app_ownership::select(tnbp::GetSiteIndexFromBond(edges),options.ownerfile,comm);
    init_random_tensor_state(ctx,edges,
		options.p_bond_dim,
		options.v_bond_dim,
		W,site_idx,site_to_mpi_rank,
		rng,comm);
  } else {
    auto filename = tnbp::make_rank_filename(options.inputname,
					     mpi_rank);
    std::ifstream ifs(filename,std::ios::binary);
    if(!ifs) throw std::runtime_error("cannot open checkpoint: "+filename);
    std::vector<Tensor> E;
    std::vector<int> edge_idx;
    tnbp::LoadTPS(ctx,ifs,W,site_idx,site_to_mpi_rank,E,edge_idx);
    if(!ifs) throw std::runtime_error("failed to read checkpoint: "+filename);
    if(mpi_rank==0)
      std::cout << " " << tnbp::make_timestamp()
                << " loaded wavefunction data" << std::endl;
  }

  if(!options.inputname.empty() && !options.ownerfile.empty()) {
    auto requested=app_ownership::select(tnbp::GetSiteIndexFromBond(edges),options.ownerfile,comm);
    app_ownership::check_legacy(site_to_mpi_rank,requested,comm);
  }
  app_ownership::print(site_to_mpi_rank,options.ownerfile.empty() ?
    (options.inputname.empty() ? "automatic" : "checkpoint") : options.ownerfile,comm);

  // Physical dimensions are local and may vary by site. An empty rank needs
  // no representative tensor and simply keeps these arrays empty.
  if(W.size()!=site_idx.size())
    throw std::invalid_argument("tensor count differs from local site labels");
  std::vector<int> physical_dims;
  physical_dims.reserve(W.size());
  for(std::size_t i=0;i<W.size();++i) {
    const auto shape=tcapi::shape(ctx,W[i]);
    const auto degree=tnbp::GetSurroundingBondIndex(site_idx[i],edges).size();
    if(shape.size()!=degree+1 || shape.back()<1 ||
       shape.back()>std::numeric_limits<int>::max())
      throw std::invalid_argument("invalid tensor physical axis at site "+std::to_string(site_idx[i]));
    physical_dims.push_back(static_cast<int>(shape.back()));
  }

  tnbp::MCMCWorkspace<Tensor,int> work;
  auto& lines = work.lines;
  std::vector<std::vector<int>> site_states(num_batches);
  for(std::size_t batch=0; batch < num_batches; batch++) {
    site_states[batch].resize(site_idx.size());
    for(std::size_t i=0;i<site_idx.size();++i) {
      std::uniform_int_distribution<int> distribution(0,physical_dims[i]-1);
      site_states[batch][i]=distribution(rng);
    }
#ifdef DEBUG
    auto sites = tnbp::GetSiteIndexFromBond(edges);
    std::vector<int> state_send(sites.size());
    std::vector<int> state_recv(sites.size());
    std::fill(state_send.begin(),state_send.end(),0);
    std::fill(state_recv.begin(),state_recv.end(),0);
    for(std::size_t adrs=0; adrs < site_idx.size(); adrs++) {
      auto it_site_adrs = std::find(sites.begin(),sites.end(),
				     site_idx[adrs]);
      auto site_adrs = std::distance(sites.begin(),it_site_adrs);
      state_send[site_adrs] = site_states[batch][adrs];
    }
    MPI_Reduce(state_send.data(),state_recv.data(),
	       static_cast<int>(sites.size()),MPI_INT,MPI_SUM,0,comm);
    if( mpi_rank == 0 ) {
      std::cout << " " << tnbp::make_timestamp()
		<< " site_states =";
      for(const auto & state : state_recv) {
	std::cout << state;
      }
      std::cout << std::endl;
    }
#endif
  }

  if( options.linefile.empty() ) {
    tnbp::bmps_setup_lines(edges,options.init_line,lines);
  } else {
    if( mpi_rank == 0 ) {
      // load lines
      lines = read_integer_table<int>(options.linefile);
    }
    bcast_lines(lines,0,comm,MPI_INT);
  }

  auto& manual_groups = work.manual_grouping;
  if(!options.groupfile.empty()) {
    // Root reads once; flatten with target labels so all ranks use identical groups.
    std::vector<std::vector<int>> rows;
    if(mpi_rank==0) {
      std::ifstream input(options.groupfile);
      if(!input) throw std::runtime_error("cannot open grouping file: "+options.groupfile);
      manual_groups=tnbp::read_mcmc_grouping<int>(input,lines);
      for(const auto& line:manual_groups) for(const auto& group:line) rows.push_back(group);
    }
    bcast_lines(rows,0,comm,MPI_INT);
    manual_groups.resize(lines.size());
    std::size_t k=0;
    for(std::size_t l=0;l<lines.size();++l) {
      manual_groups[l].resize(lines[l].size());
      for(auto& group:manual_groups[l]) group=rows.at(k++);
    }
    tnbp::validate_mcmc_grouping(edges,lines,manual_groups);
  }

  if( mpi_rank == 0 ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " end initialization of random site_states " << std::endl;
  }

  // Initialize all batches with one shared grouping and workspace.
  tnbp::mcmc_initialize(ctx,edges,W,site_idx,site_to_mpi_rank,
                       site_states,options.bmps_max_bp_steps,
                       options.bmps_bp_tolerance,options.bmps_bond_dim,
                       options.bmps_sv_min,options.bmps_tg_err,work,comm);

  if( mpi_rank == 0 ) {
    std::cout << " " << tnbp::make_timestamp()
	      << " mcmc_initialization is finished " << std::endl;
  }

  std::size_t min_line_idx = 0;
  std::size_t max_line_idx = (lines.size()>0) ? (lines.size()-1) : 0;
  if( mpi_rank == 0 ) {
    for(std::size_t line_idx=min_line_idx; line_idx <= max_line_idx; line_idx++) {
      std::cout << " " << tnbp::make_timestamp ()
		<< " site in line " << line_idx << ":";
      for(std::size_t site_adrs=0; site_adrs < lines[line_idx].size(); site_adrs++) {
	std::cout << " " << lines[line_idx][site_adrs];
      }
      std::cout << std::endl;
    }
  }

#ifdef DEBUG
  for(std::size_t batch=0; batch < num_batches; batch++) {
    for(int rank=0; rank < mpi_size; rank++) {
      if( mpi_rank == rank ) {
	for(std::size_t adrs=0; adrs < site_idx.size(); adrs++) {
	  auto shape_v = tcapi::shape(ctx,work.V[batch][adrs]);
	  std::cout << " " << tnbp::make_timestamp()
		    << " tensor at site " << site_idx[adrs] << ":";
	  for(const auto & d : shape_v) {
	    std::cout << " " << d;
	  }
	  std::cout << std::endl;
	}
      }
      MPI_Barrier(comm);
    }
  }
#endif

  std::vector<std::vector<std::vector<int>>> sampling_states(
      options.num_samples,
      std::vector<std::vector<int>>(num_batches,std::vector<int>(site_idx.size())));
  if(!sampling_states.empty()) sampling_states.front() = site_states;

  tnbp::mcmc_sampling(ctx,edges,W,site_idx,site_to_mpi_rank,
                      sampling_states,options.sampling_method,rng,
                      options.bmps_max_bp_steps,options.bmps_bp_tolerance,
                      options.bmps_bond_dim,options.bmps_sv_min,options.bmps_tg_err,
                      work,comm);

  if( mpi_rank == 0 ) {
    std::cout << " # sample sample_id.batch : bitstring " << std::endl;
  }
  auto sites = tnbp::GetSiteIndexFromBond(edges);
  std::vector<int> state_at_site_send(sites.size());
  std::vector<int> state_at_site_recv(sites.size());
  for(std::size_t sample=0; sample < static_cast<std::size_t>(options.num_samples); sample++) {
    for(std::size_t batch=0; batch < static_cast<std::size_t>(num_batches); batch++) {
      std::fill(state_at_site_send.begin(),state_at_site_send.end(),0);
      std::fill(state_at_site_recv.begin(),state_at_site_recv.end(),0);
      for(std::size_t adrs=0; adrs < site_idx.size(); adrs++) {
	auto it_site_adrs = std::find(sites.begin(),sites.end(),
				      site_idx[adrs]);
	auto site_adrs = std::distance(sites.begin(),it_site_adrs);
	state_at_site_send[site_adrs] = sampling_states[sample][batch][adrs];
      }
      MPI_Reduce(state_at_site_send.data(),state_at_site_recv.data(),
		 static_cast<int>(sites.size()),MPI_INT,MPI_SUM,
		 0,comm);
      if( mpi_rank == 0 ) {
	std::cout << " sample " << sample << "." << batch << ": ";
	for(const auto & state : state_at_site_recv) {
	  std::cout << state;
	}
	std::cout << std::endl;
      }
    }
  }

  } catch(const std::exception& error) {
    std::cerr << "sampling: " << error.what() << std::endl;
    MPI_Abort(comm,1);
  }
  MPI_Finalize();
  return 0;

}
