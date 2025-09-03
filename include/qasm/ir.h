// qasm/ir.hpp
#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <cstddef>

namespace qasm {

struct QReg { std::string name; std::size_t size{}; };
struct CReg { std::string name; std::size_t size{}; };

enum class Op {
  U3,U2,U1,RX,RY,RZ,H,X,Y,Z,S,SDG,T,TDG,ID,CX,CZ,SWAP,CCX,CSWAP,
  MEASURE,BARRIER,RESET,CUSTOM
};

struct QArg { std::string reg; std::size_t index{}; };

struct Instruction {
  Op op{};
  std::string name;               // for CUSTOM 
  std::vector<double> params;
  std::vector<QArg> qubits;
  bool has_classical{false};      // for measure 
  std::string c_reg; std::size_t c_index{};
  std::size_t line{};
  std::string raw;
};

struct Program {
  std::vector<QReg> qregs;
  std::vector<CReg> cregs;
  std::vector<Instruction> instructions;
};

struct ParseError : std::runtime_error {
  std::size_t line, col;
  ParseError(std::size_t l, std::size_t c, const std::string& msg)
  : std::runtime_error("QASM parse error ("+std::to_string(l)+":"+std::to_string(c)+"): "+msg),
    line(l), col(c) {}
};

} // namespace qasm
