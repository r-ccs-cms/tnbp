#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Option {
  std::string backend("ibm_kobe");
  std::string circuit("circuit.qasm");
  std::vector<int> num_gates;
  size_t max_bp_iterations = 50;
};

Option generate_options(int argc, char *argv[]) {

  Option option;
  for(int i=0; i < argc; i++) {
    if ( std::string(argv[i]) == "--backend" ) {
      option.backend = std::string(argv[++i]);
    }
    if ( std::string(argv[i]) == "--circuit" ) {
      option.circuit = std::string(argv[++i]);
    }
    if ( std::string(argv[i]) == "--num_gates" ) {
      std::stringstream ss(argv[++i]);
      std::string token;
      while (std::getline(ss, token, ',')) {
	option.num_gates.push_back(std::stoi(token));
      }
    }
    if ( std::string(argv[i]) == "--max_bp_iterations" ) {
      option.max_bp_iterations = std::atoi(argv[++i]);
    }
  }
  return option;
}
