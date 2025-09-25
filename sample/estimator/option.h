#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Option {

  // inputs for problem settings
  std::string backend = "ibm_kobe";
  std::string circuit = "circuit.qasm";
  std::string sparsepauli = "sparsepauliop.txt";

  // inputs for controlled parameters
  std::vector<int> num_gates;
  size_t max_bp_iterations = 50;
  size_t max_bond_dim = 100;
  double bp_tolerance = 1.0e-4;
  double sv_min = 1.0e-8;
  double truncation_error = 1.0e-4;
  
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
    if ( std::string(argv[i]) == "--bp_tolerance" ) {
      option.bp_tolerance = std::atof(argv[++i]);
    }
    if ( std::string(argv[i]) == "--max_bond_dim" ) {
      option.max_bond_dim = std::atoi(argv[++i]);
    }
    if ( std::string(argv[i]) == "--sv_min" ) {
      option.sv_min = std::atof(argv[++i]);
    }
    if ( std::string(argv[i]) == "--truncation_error" ) {
      option.truncation_error = std::atof(argv[++i]);
    }
    if ( std::string(argv[i]) == "--sparse_pauli" ) {
      option.sparsepauli = std::string(argv[++i]);
    }
  }
  return option;
}

void cout_options(const Option & option) {
  std::cout.precision(16);
  std::cout << "# backend: " << option.backend << std::endl;
  std::cout << "# circuit: " << option.circuit << std::endl;
  std::cout << "# sparsepauli: " << option.sparsepauli << std::endl;
  std::cout << "# num gates:";
  for(const auto ng : option.num_gates) {
    std::cout << " " << ng;
  }
  std::cout << std::endl;
  std::cout << "# max_bp_iterations: " << option.max_bp_iterations << std::endl;
  std::cout << "# bp_tolerance: " << option.bp_tolerance << std::endl;
  std::cout << "# max_bond_dim: " << option.max_bond_dim << std::endl;
  std::cout << "# sv_min: " << option.sv_min << std::endl;
  std::cout << "# truncation_error: " << option.truncation_error << std::endl;
}
