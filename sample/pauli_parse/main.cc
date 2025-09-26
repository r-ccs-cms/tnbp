// main.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <cassert>
#include "pauli/sparse_pauli.h"
#include "pauli/pauli_string.h"

#ifdef USE_REAL
using ElemT = double;
#else
using ElemT = std::complex<double>;
#endif

static void print_localops_compact(const std::vector<pauli::LocalOp> & ops) {
  std::cout << "(";
  for(size_t i = 0; i < ops.size(); ++i) {
    std::cout << "q" << ops[i].qubit << ":" << ops[i].op;
    if( i+1 < ops.size() ) std::cout << " ";
  }
  std::cout << ")";
}

int main(int argc, char** argv){
  if ( argc < 2 ) {
    std::cerr << "Usage: " << argv[0] << " <sparsepauli_file.txt> [--show-all]\n"
	      << "  <sparse_pauli_op.txt> : lines of 'real imag PauliString' or 'real PauliString'\n"
	      << "  --show-all            : (optional) print Identity-included LocalOps (can be huge)\n";
    return 1;
  }

  const std::string filename = argv[1];
  const bool show_all = (argc >= 3 && std::string(argv[2]) == "--show-all");

  try {
    auto terms = pauli::load_sparse_pauli_op<ElemT>(filename);
    std::cout << "Loaded " << terms.size() << " Pauli terms from " << filename << std::endl;

    assert(!terms.empty());
    std::size_t n_qubits = terms.front().pauli_string.size();

    const std::size_t MAX_SHOW = 10;
    std::size_t shown = 0;
    
    for(const auto & t : terms) {

      if ( shown >= MAX_SHOW ) {
	if ( shown >= MAX_SHOW ) {
	  std::cout << "... (truncated; showing first " << MAX_SHOW << " terms)\n";
	  break;
	}
      }
      
      std::cout << "coeff = " << t.coeff
		<< ", pauli_string = " << t.pauli_string << std::endl;

      const auto & c = t.coeff;
      const auto & s = t.pauli_string;

      auto all_ops = pauli::string_to_localops(s);
      auto nz_ops  = pauli::string_to_non_identity_localops(s);

      assert(all_ops.size() == n_qubits);

      std::cout << "term[" << shown << "]: "
		<< std::setprecision(12)
		<< "(" << c << ") * " << s << std::endl;

      std::cout << "  non-identity localops: ";
      print_localops_compact(nz_ops);
      std::cout << std::endl;

      if( show_all ) {
	std::cout << "  all localops (|ops|=" << all_ops.size() << "): ";
	print_localops_compact(all_ops);
	std::cout << std::endl;
      }
    }
  } catch (const std::exception & ex) {
    std::cerr << "Error reading file: " << ex.what() << std::endl;
    return 1;
  }
  
  return 0;
}
