#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <complex>
#include <cctype>

namespace pauli {

  // Term
  template <class CoeffT>
  struct Term {
    CoeffT coeff;
    std::string pauli_string; // e.g. "ZIIX"
  };

  // delete spaces 
  inline void trim_inplace(std::string& s) {
    auto is_ns = [](unsigned char c){ return !std::isspace(c); };
    auto b = std::find_if(s.begin(), s.end(), is_ns);
    auto e = std::find_if(s.rbegin(), s.rend(), is_ns).base();
    if (b < e) s.assign(b, e); else s.clear();
  }
  
  // get extension (no dot / e.g. "txt")
  inline std::string get_extension(const std::string& filename) {
    auto pos = filename.find_last_of('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos + 1);
  }
  
  // Generate CoeffT (rea/imag -> real or complex)
  template <class CoeffT>
  inline CoeffT make_coeff(double re, double im) {
    if constexpr (std::is_same_v<CoeffT, std::complex<double>>) {
      return CoeffT(re, im);
    } else {
      // if real, neglect imaginary part
      return static_cast<CoeffT>(re);
    }  
  }
  
  // Load text compatible to SparsePauliOp of Qiskit
  // Format:
  //   - "real imag PauliString"  e.g.: 0.0 1.2 ZZI
  //   - "real PauliString"       e.g.: 0.5 ZII
  // Comment: Ignore the lines starting with '#' or '//'. Also ignore the empty line.
  template <class CoeffT = std::complex<double>>
  inline std::vector<Term<CoeffT>> load_sparse_pauli_op(const std::string& filename) {
    std::string ext = get_extension(filename);
    if (ext != "txt") {
      throw std::runtime_error("Unsupported file extension: " + ext + " (only .txt is supported for now)");
    }
    
    std::ifstream in(filename);
    if (!in) {
      throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::vector<Term<CoeffT>> op;
    std::string line;
    size_t lineno = 0;
    
    while (std::getline(in, line)) {
        ++lineno;
        trim_inplace(line);
        if (line.empty()) continue;
        if (line.rfind("#", 0) == 0) continue;       // beginning line '#'
        if (line.rfind("//", 0) == 0) continue;      // beginning line '//'
	
        std::istringstream iss(line);
        double re = 0.0, im = 0.0;
        std::string pauli;
	
        // try "re im pauli"
        if ( (iss >> re >> im >> pauli) ) {
	  op.push_back({ make_coeff<CoeffT>(re, im), pauli });
	  continue;
        }
	
        // if it is not possible, try "re pauli"
        iss.clear();
        iss.str(line);
        if ( (iss >> re >> pauli) ) {
	  op.push_back({ make_coeff<CoeffT>(re, 0.0), pauli });
	  continue;
        }
	
        throw std::runtime_error("Parse error at line " + std::to_string(lineno) + ": \"" + line + "\"");
    }
    
    return op;
  }
  
} // namespace pauli
