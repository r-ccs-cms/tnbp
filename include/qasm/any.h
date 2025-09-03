// qasm/any.hpp
#pragma once
#include "v2/parser.hpp"
#include "v3/parser.hpp"

namespace qasm {

inline Program parse_any(std::string_view src) {
  // Simple header sniff
  auto pos = src.find("OPENQASM");
  if (pos != std::string_view::npos) {
    auto tail = src.substr(pos, 32);
    if (tail.find('3') != std::string_view::npos) {
      return v3::Parser{}.parse(src);
    }
    // default to v2 if header says 2 or unknown
    return v2::Parser{}.parse(src);
  }
  // No header? assume v2 for compatibility
  return v2::Parser{}.parse(src);
}

} // namespace qasm
