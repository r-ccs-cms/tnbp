#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Option {

  // inputs for problem settings
  std::string inputname;
  std::string outputname;

  bool do_contraction = false;

  std::size_t state_bond_dim = 1;
  uint32_t seed = 1729;

  // boundary-mps parameters
  std::size_t bmps_bond_dim = 1;
  Real bmps_bp_tolerance = 1.0e-8;
  int bmps_max_bp_steps = 10;
  Real bmps_tg_err = 1.0e-8;
  Real bmps_sv_min = 1.0e-8;

  // initial line to define the lines
  std::vector<int> init_line;

};

Option generate_options(int argc, char *argv[]) {

  Option option;
  for(int i=0; i < argc; i++) {
    if ( std::string(argv[i]) == "--inputname" ) {
      option.inputname = std::string(argv[++i]);
    }
    if ( std::string(argv[i]) == "--outputname" ) {
      option.outputname = std::string(argv[++i]);
    }
    if ( std::string(argv[i]) == "--do_contraction" ) {
      if( std::atoi(argv[++i]) != 0 ) {
	option.do_contraction = true;
      }
    }
    if( std::string(argv[i]) == "--state_bond_dim" ) {
      option.state_bond_dim = static_cast<size_t>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--seed" ) {
      option.seed = static_cast<uint32_t>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--bmps_bond_dim" ) {
      option.bmps_bond_dim = static_cast<size_t>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--bmps_bp_tolerance" ) {
      option.bmps_bp_tolerance = static_cast<Real>(std::atof(argv[++i]));
    }
    if( std::string(argv[i]) == "--bmps_max_bp_steps" ) {
      option.bmps_max_bp_steps = static_cast<int>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--bmps_tg_err" ) {
      option.bmps_tg_err = static_cast<Real>(std::atof(argv[++i]));
    }
    if( std::string(argv[i]) == "--bmps_sv_min" ) {
      option.bmps_sv_min = static_cast<Real>(std::atof(argv[++i]));
    }
    if( std::string(argv[i]) == "--init_line" ) {
      std::stringstream ss(argv[++i]);
      std::string token;
      while (std::getline(ss,token,',')) {
	option.init_line.push_back(static_cast<int>(std::stoi(token)));
      }
    }
  }
  return option;
}

void cout_options(const Option & option) {
  std::cout.precision(16);
  std::cout << "# input graph file name: " << option.inputname << std::endl;
  if( option.do_contraction ) {
    std::cout << "# do_contraction: true" << std::endl;
    std::cout << "# bond dimension: " << option.state_bond_dim << std::endl;
    std::cout << "# seed for random tensor contraction: " << option.seed << std::endl;
  }
  std::cout << "# bond dimsion for boundary-mps: " << option.bmps_bond_dim << std::endl;
  std::cout << "# tolerance for belief propagation: " << option.bmps_bp_tolerance << std::endl;
  std::cout << "# max belief propagation steps: " << option.bmps_max_bp_steps << std::endl;
  std::cout << "# target error of boundary-mps truncation: " << option.bmps_tg_err << std::endl;
  std::cout << "# minimum of singular value in boundary-mps truncation: " << option.bmps_sv_min << std::endl;
  std::cout << "# output file name: " << option.outputname << std::endl;
  if( !option.init_line.empty() ) {
    std::cout << "# init line:";
    for(auto const & site : option.init_line) {
      std::cout << " " << site;
    }
    std::cout << std::endl;
  }
  
}
