#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Option {
  int num_qubits = 4;

  double Jz = 1.0;
  double hz = 0.0;
  double hx = 1.0;
  double dt = 0.01 * PI_v<double>;

  int num_steps = 1;
  
  size_t max_bp_iterations = 50;
  size_t max_bond_dim = 100;
  double bp_tolerance = 1.0e-4;
  double sv_min = 1.0e-8;
  double truncation_error = 1.0e-4;

};

Option generate_options(int argc, char *argv[]) {

  Option option;
  for(int i=0; i < argc; i++) {
    if ( std::string(argv[i]) == "--num_qubits" ) {
      option.num_qubits = std::atoi(argv[++i]);
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
    if ( std::string(argv[i]) == "--Jz" ) {
      option.Jz = std::atof(argv[++i]);
    }
    if ( std::string(argv[i]) == "--hz" ) {
      option.hz = std::atof(argv[++i]);
    }
    if ( std::string(argv[i]) == "--hx" ) {
      option.hx = std::atof(argv[++i]);
    }
    if ( std::string(argv[i]) == "--dt" ) {
      option.dt = std::atof(argv[++i]) * PI_v<double>;
    }
    if ( std::string(argv[i]) == "--num_steps" ) {
      option.num_steps = std::atoi(argv[++i]);
    }
  }
  return option;
}

void cout_options(const Option & option) {
  std::cout.precision(16);
  std::cout << "# num_qubits: " << option.num_qubits << std::endl;
  std::cout << "# Jz: " << option.Jz << std::endl;
  std::cout << "# hz: " << option.hz << std::endl;
  std::cout << "# hx: " << option.hx << std::endl;
  std::cout << "# dt: " << option.dt << std::endl;
  std::cout << "# num_steps: " << option.num_steps << std::endl;
  std::cout << "# max_bp_iterations: " << option.max_bp_iterations << std::endl;
  std::cout << "# bp_tolerance: " << option.bp_tolerance << std::endl;
  std::cout << "# max_bond_dim: " << option.max_bond_dim << std::endl;
  std::cout << "# sv_min: " << option.sv_min << std::endl;
  std::cout << "# truncation_error: " << option.truncation_error << std::endl;
}
