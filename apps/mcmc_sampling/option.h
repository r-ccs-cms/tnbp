#pragma once
#include <cstddef>
#include <cstdint>
#include <charconv>
#include <cmath>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Option {
  bool help = false;

  // inputs for problem settings
  std::string ownerfile;
  std::string portable_input;
  std::string inputname;
  std::string inputedge;
  std::string outputname;
  std::uint32_t seed = 1729;

  // initialization for random initial state
  std::size_t p_bond_dim = 1;
  std::size_t v_bond_dim = 1;

  // boundary-mps parameters
  std::size_t bmps_bond_dim = 1;
  Real bmps_bp_tolerance = 1.0e-8;
  int bmps_max_bp_steps = 10;
  Real bmps_tg_err = 1.0e-8;
  Real bmps_sv_min = 1.0e-8;

  // initial line to define the lines
  std::vector<int> init_line;
  std::string linefile;
  std::string groupfile;

  // parameters for sampling
  int num_samples = 1;
  int num_batches = 1;
  int sampling_method = 2;
  int update_sites = 1;

};

/** @brief Print usage, defaults and sampling restrictions without loading inputs. */
inline void cout_help() {
  std::cout << R"HELP(Usage: mcmc_sampling[_complex] --inputedge FILE (--init_line IDS | --linefile FILE) [options]

Options (defaults in brackets):
  -h, --help                   [—]
      Show this help and exit successfully; no input files are required.
  --inputedge FILE             [required]
      Undirected graph edge file, one pair of site IDs per row.
  --ownerfile FILE             [none]
      Ownership rows: site_id mpi_rank; omitted uses the automatic partition.
  --portable_input DIRECTORY   [none]
      Load a portable TPS on a new MPI partition, including its graph.
  --inputname PREFIX           [none]
      Load rank-specific TPS files; otherwise generate random tensors.
  --outputname NAME            [none]
      Legacy option: displayed only; does not write a file. Redirect stdout.
  --seed N                     [1729]
      Random seed for tensor generation, initial states and sampling.
  --p_bond_dim N               [1]
      Physical dimension of random tensors; 1 gives fixed state zero.
  --v_bond_dim N               [1]
      Virtual bond dimension of random tensors.
  --init_line IDS              [none]
      Comma-separated initial boundary IDs for automatic line construction.
  --linefile FILE              [none]
      Explicit lines: one row of space-separated site IDs per line. Overrides --init_line.
  --groupfile FILE             [none]
      Manual three-line groups; omitted means automatic grouping on the original graph.
  --bmps_bond_dim N            [1]
      Maximum retained boundary-MPS bond dimension.
  --bmps_bp_tolerance X        [1e-8]
      Belief-propagation residual stopping tolerance.
  --bmps_max_bp_steps N        [10]
      Maximum BP iterations per boundary truncation.
  --bmps_tg_err X              [1e-8]
      Target boundary truncation error.
  --bmps_sv_min X              [1e-8]
      Singular-value cutoff for boundary truncation.
  --num_samples N              [1]
      Saved samples per batch; one sample follows a full forward/backward sweep.
  --num_batches N              [1]
      Number of chains initialized and sampled.
  --sampling_method N          [2]
      0: Metropolis; 1: pairwise Barker; 2: full-state roulette heat bath.
  --update_sites N             [1]
      1 or 2 consecutive line sites; 2 requires method 2. Singleton lines are skipped.

Both line modes require at least two nonempty lines covering all graph sites.
--inputname ignores the random-tensor dimensions; checkpoint backend, scalar type
and MPI partition must match. Preserve the saved incident-edge order in FILE;
duplicate edges are rejected. Empty ranks and per-site physical dimensions are supported.
Dimensions/batches must be positive; samples/BP steps may be zero. Tolerances
must be finite and nonnegative; malformed values and unknown options are rejected.
Samples are printed to stdout in site-label order.
No burn-in or thinning is performed automatically. For uniform binary states,
use the default method 2 to avoid Metropolis round-trip periodicity.
Method 2 enumerates every local state (all d_i*d_j pairs for two sites), without
charge filtering. The current state is included in the roulette draw.

Example:
  mpirun -np 1 mcmc_sampling --inputedge chain.txt --init_line 0 \
    --p_bond_dim 2 --v_bond_dim 1 --num_samples 100
)HELP";
}

Option generate_options(int argc, char *argv[]) {

  Option option;
  // Handle help before parsing values or touching input files.
  for(int i=1; i<argc; ++i) {
    const std::string arg(argv[i]);
    if(arg=="--help" || arg=="-h") { option.help=true; return option; }
  }
  auto integer=[](const std::string& flag, const std::string& value,
                  long long minimum, long long maximum) {
    long long parsed=0;
    const auto result=std::from_chars(value.data(),value.data()+value.size(),parsed);
    if(result.ec!=std::errc{} || result.ptr!=value.data()+value.size() ||
       parsed<minimum || parsed>maximum)
      throw std::invalid_argument(flag+" requires an integer in ["+
          std::to_string(minimum)+", "+std::to_string(maximum)+"]");
    return parsed;
  };
  auto real=[](const std::string& flag, const std::string& value) {
    std::size_t used=0;
    double parsed;
    try { parsed=std::stod(value,&used); }
    catch(const std::exception&) {
      throw std::invalid_argument(flag+" requires a finite nonnegative number");
    }
    if(used!=value.size() || !std::isfinite(parsed) || parsed<0 ||
       parsed>std::numeric_limits<Real>::max())
      throw std::invalid_argument(flag+" requires a finite nonnegative number");
    return static_cast<Real>(parsed);
  };
  const auto int_max=std::numeric_limits<int>::max();
  for(int i=1; i<argc; ++i) {
    const std::string flag(argv[i]);
    if(i+1>=argc || std::string(argv[i+1]).rfind("--",0)==0 ||
       std::string(argv[i+1])=="-h")
      throw std::invalid_argument(flag+" requires a value");
    const std::string value(argv[++i]);
    if(value.empty()) throw std::invalid_argument(flag+" requires a nonempty value");
    if(flag=="--ownerfile") option.ownerfile=value;
    else if(flag=="--portable_input") option.portable_input=value;
    else if(flag=="--inputname") option.inputname=value;
    else if(flag=="--inputedge") option.inputedge=value;
    else if(flag=="--outputname") option.outputname=value;
    else if(flag=="--linefile") option.linefile=value;
    else if(flag=="--groupfile") option.groupfile=value;
    else if(flag=="--seed") option.seed=static_cast<std::uint32_t>(
        integer(flag,value,0,std::numeric_limits<std::uint32_t>::max()));
    else if(flag=="--p_bond_dim") option.p_bond_dim=integer(flag,value,1,int_max);
    else if(flag=="--v_bond_dim") option.v_bond_dim=integer(flag,value,1,int_max);
    else if(flag=="--bmps_bond_dim") option.bmps_bond_dim=integer(flag,value,1,int_max);
    else if(flag=="--bmps_bp_tolerance") option.bmps_bp_tolerance=real(flag,value);
    else if(flag=="--bmps_tg_err") option.bmps_tg_err=real(flag,value);
    else if(flag=="--bmps_sv_min") option.bmps_sv_min=real(flag,value);
    else if(flag=="--bmps_max_bp_steps") option.bmps_max_bp_steps=integer(flag,value,0,int_max);
    else if(flag=="--num_samples") option.num_samples=integer(flag,value,0,int_max);
    else if(flag=="--num_batches") option.num_batches=integer(flag,value,1,int_max);
    else if(flag=="--sampling_method") option.sampling_method=integer(flag,value,0,2);
    else if(flag=="--update_sites") option.update_sites=integer(flag,value,1,2);
    else if(flag=="--init_line") {
      if(value.back()==',') throw std::invalid_argument(flag+" contains an empty site ID");
      std::stringstream stream(value);
      std::string token;
      while(std::getline(stream,token,','))
        option.init_line.push_back(static_cast<int>(integer(flag,token,
            std::numeric_limits<int>::min(),int_max)));
    } else throw std::invalid_argument("unknown option: "+flag);
  }
  if(option.update_sites==2 && option.sampling_method!=2)
    throw std::invalid_argument("--update_sites 2 requires --sampling_method 2");
  return option;
}

void cout_options(const Option & option) {
  std::cout.precision(16);
  if( option.inputname.empty() ) {
    std::cout << "# initial tensor network is random initial state"
	      << std::endl;
    std::cout << "# physical bond dimension of random inital state: "
	      << option.p_bond_dim << std::endl;
    std::cout << "# physical bond dimension of random inital state: "
	      << option.v_bond_dim << std::endl;
  } else {
    std::cout << "# input tensor network data name: "
	      << option.inputname << std::endl;
  }
  std::cout << "# seed for random tensor contraction: "
	    << option.seed << std::endl;
  std::cout << "# bond dimsion for boundary-mps: "
	    << option.bmps_bond_dim << std::endl;
  std::cout << "# tolerance for belief propagation: "
	    << option.bmps_bp_tolerance << std::endl;
  std::cout << "# max belief propagation steps: "
	    << option.bmps_max_bp_steps << std::endl;
  std::cout << "# target error of boundary-mps truncation: "
	    << option.bmps_tg_err << std::endl;
  std::cout << "# minimum of singular value in boundary-mps truncation: "
	    << option.bmps_sv_min << std::endl;
  std::cout << "# output file name: "
	    << option.outputname << std::endl;
  if( !option.init_line.empty() ) {
    std::cout << "# init line:";
    for(auto const & site : option.init_line) {
      std::cout << " " << site;
    }
    std::cout << std::endl;
  }
  if( !option.linefile.empty() ) {
    std::cout << "# file for lines: " << option.linefile << std::endl;
  }
  if(!option.groupfile.empty()) std::cout << "# grouping file: " << option.groupfile << std::endl;
  std::cout << "# number of samples: " << option.num_samples << std::endl;
  std::cout << "# number of batches: " << option.num_batches << std::endl;
  std::cout << "# update sites: " << option.update_sites << std::endl;
  std::cout << "# sampling method: " << option.sampling_method << std::endl;
#ifdef USE_COMPLEX
  std::cout << "# complex version " << std::endl;
#endif

}
