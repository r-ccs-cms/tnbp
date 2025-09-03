// main.cpp
#include <iostream>
#include "qasm/any.hpp"

static void dump(const qasm::Program& p){
  std::cout << "qregs:\n";
  for (auto& r: p.qregs) std::cout << "  qreg " << r.name << "[" << r.size << "]\n";
  std::cout << "cregs:\n";
  for (auto& r: p.cregs) std::cout << "  creg " << r.name << "[" << r.size << "]\n";
  std::cout << "instructions: " << p.instructions.size() << "\n";
  for (auto& ins: p.instructions) {
    std::cout << "  line " << ins.line << ": ";
    if (ins.op==qasm::Op::CUSTOM) std::cout << ins.name;
    else std::cout << (int)ins.op;
    if (!ins.params.empty()) {
      std::cout << " (";
      for (size_t i=0;i<ins.params.size();++i){ if(i) std::cout<<", "; std::cout<<ins.params[i]; }
      std::cout << ")";
    }
    std::cout << "  ";
    for (size_t i=0;i<ins.qubits.size();++i){
      if (i) std::cout<<", ";
      std::cout << ins.qubits[i].reg << "[" << ins.qubits[i].index << "]";
    }
    if (ins.op==qasm::Op::MEASURE) std::cout << " -> " << ins.c_reg << "[" << ins.c_index << "]";
    std::cout << "\n";
  }
}

int main(){
  const char* qasm2 = R"(OPENQASM 2.0;
include "qelib1.inc";
qreg q[3]; creg c[3];
h q;
cx q[0], q[1];
rz(pi/2) q[2];
measure q[1] -> c[1];
barrier q[0], q[2];
)";

  const char* qasm3 = R"(OPENQASM 3;
include "stdgates.inc";
qubit[4] q; bit[4] c;
h q[0:3];
cx q[0], q[3];
c[1] = measure q[1];
bit x = measure q[2];
barrier q;
reset q[0];
)";

  try{
    std::cout << "=== V2 ===\n";
    dump(qasm::parse_any(qasm2));
    std::cout << "\n=== V3 ===\n";
    dump(qasm::parse_any(qasm3));
  } catch(const qasm::ParseError& e){
    std::cerr << e.what() << "\n";
    return 1;
  }
}
