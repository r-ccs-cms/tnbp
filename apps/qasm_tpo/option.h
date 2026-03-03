#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Option {

  // inputs for problem settings
  std::string circuit = "circuit.qasm";

  // cutoff function
  bool do_opt_tpo = false;
  double opt_tpo_eps = 1.0e-8;
  

};

Option generate_options(int argc, char *argv[]) {

  Option option;
  for(int i=0; i < argc; i++) {
    if ( std::string(argv[i]) == "--circuit" ) {
      option.circuit = std::string(argv[++i]);
    }
    if ( std::string(argv[i]) == "--do_opt_tpo" ) {
      if( std::atoi(argv[++i]) != 0 ) {
	option.do_opt_tpo = true;
      }
    }
    if ( std::string(argv[i]) == "--opt_tpo_eps" ) {
      option.opt_tpo_eps = std::atof(argv[++i]);
    }
  }
  return option;
}

void cout_options(const Option & option) {
  std::cout.precision(16);
  std::cout << "# circuit: " << option.circuit << std::endl;
  if( option.do_opt_tpo ) {
    std::cout << "# do_opt_tpo: true" << std::endl;
    std::cout << "# opt_tpo_eps: " << option.opt_tpo_eps << std::endl;
  }
}
