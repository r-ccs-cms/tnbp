#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Option {

  // inputs for problem settings
  std::string inputname;
  std::string outputname;

  //
  bool do_contraction = false;

  std::size_t bond_dim = 1;
  std::size_t num_fin_edges = 3;
  std::string schedule_type = "random";
  uint32_t seed = 1729;
  uint32_t seed_for_ten = 1729;
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
    if( std::string(argv[i]) == "--bond_dim" ) {
      option.bond_dim = static_cast<size_t>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--schedule_type" ) {
      option.schedule_type = std::string(argv[++i]);
    }
    if( std::string(argv[i]) == "--seed" ) {
      option.seed = static_cast<uint32_t>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--num_fin_edges" ) {
      option.num_fin_edges = static_cast<std::size_t>(std::atoi(argv[++i]));
    }
    if( std::string(argv[i]) == "--seed_for_ten" ) {
      option.seed_for_ten = static_cast<uint32_t>(std::atoi(argv[++i]));
    }
  }
  return option;
}

void cout_options(const Option & option) {
  std::cout.precision(16);
  std::cout << "# input graph file name: " << option.inputname << std::endl;
  if( option.do_contraction ) {
    std::cout << "# do_contraction: true" << std::endl;
    std::cout << "# bond dimension: " << option.bond_dim << std::endl;
    std::cout << "# seed for random tensor contraction: " << option.seed_for_ten << std::endl;
  }
  std::cout << "# schedule type: " << option.schedule_type << std::endl;
  if( option.schedule_type == std::string("random") ) {
    std::cout << "# seed: " << option.seed << std::endl;
  }
  std::cout << "# number of final edges: " << option.num_fin_edges << std::endl;
  std::cout << "# output file name: " << option.outputname << std::endl;
}
