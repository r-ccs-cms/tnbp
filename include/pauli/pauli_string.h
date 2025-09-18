#pragma once
#include <string>
#include <vector>
#include <cctype>
#include <stdexcept>
#include <algorithm>

namespace pauli {
  
  // One-qubit local unitary operator (including identity)
  struct LocalOp {
    std::size_t qubit; // 0-based. right edge is qubit 0
    char op;           // 'I' | 'X' | 'Y' | 'Z'
  };
  
  // Example whether string is composed of I/X/Y/Z, and make all of letters capital
  inline void normalize_pauli_label(std::string& s) {
    trim_inplace(s);
    for (char& c : s) {
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      if (c != 'I' && c != 'X' && c != 'Y' && c != 'Z') {
	throw std::runtime_error("Invalid Pauli label: only I/X/Y/Z are allowed.");
      }
    }
  }
  
  // return all-site operator assuming that right-edge is qubit 0.
  // e.g.: "XIZZY" -> [{q0:'Y'}, {q1:'Z'}, {q2:'Z'}, {q3:'I'}, {q4:'X'}]
  inline std::vector<LocalOp> string_to_localops(std::string pauli_label) {
    normalize_pauli_label(pauli_label);
    const std::size_t n = pauli_label.size();
    std::vector<LocalOp> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      char u = pauli_label[n - 1 - i]; // right-edge -> left-edge
      out.push_back(LocalOp{ i, u });
    }
    return out; // qubit 0,1,2,... (ascending order)
  }
  
  // return operator composed of non-identity local operators
  // e.g.: "XIZZY" -> [{q0:'Y'}, {q1:'Z'}, {q2:'Z'}, {q4:'X'}]
  inline std::vector<LocalOp> string_to_non_identity_localops(std::string pauli_label) {
    auto all = string_to_localops(std::move(pauli_label));
    std::vector<LocalOp> out;
    out.reserve(all.size());
    for (const auto& t : all) {
      if (t.op != 'I') out.push_back(t);
    }
    return out;
  }
  
} // namespace pauli
