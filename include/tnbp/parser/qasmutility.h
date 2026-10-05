// Shared QASM instruction metadata; no tensor backend is required.
#ifndef TNBP_PARSER_QASMUTILITY_H
#define TNBP_PARSER_QASMUTILITY_H

#include "qasm/ir.h"

namespace tnbp {
  /**
     Function to get the size of qubits of the gate defined by qasm::Instruction
     @param[in] ins: qasm::Instruction corresponding to the gate
   */
  inline int OpQubitCount(const qasm::Instruction & ins) {
    switch (ins.op) {
    case qasm::Op::U3: case qasm::Op::U2: case qasm::Op::U1:
    case qasm::Op::RX: case qasm::Op::RY: case qasm::Op::RZ:
    case qasm::Op::H:  case qasm::Op::X:  case qasm::Op::Y:  case qasm::Op::Z:
    case qasm::Op::S:  case qasm::Op::SDG: case qasm::Op::T: case qasm::Op::TDG:
    case qasm::Op::ID:
      return 1;
    case qasm::Op::CX: case qasm::Op::CZ: case qasm::Op::SWAP: case qasm::Op::RZZ:
      return 2;
    case qasm::Op::CCX: case qasm::Op::CSWAP:
      return 3;
    case qasm::Op::MEASURE:
      return 1;
    case qasm::Op::RESET:
      return 1;
    case qasm::Op::BARRIER:
    case qasm::Op::CUSTOM:
      return ins.qubits.size(); // 可変長は実データ依存
    }
    return ins.qubits.size();
  }

} // namespace tnbp

#endif
