// qasm/utility.h
#pragma once

namespace qasm {

  static void dump(const Program & p){
    std::cout << "qregs:\n";
    for (auto& r: p.qregs) std::cout << "  qreg " << r.name << "[" << r.size << "]\n";
    std::cout << "cregs:\n";
    for (auto& r: p.cregs) std::cout << "  creg " << r.name << "[" << r.size << "]\n";
    std::cout << "instructions: " << p.instructions.size() << "\n";
    for (auto& ins: p.instructions) {
      std::cout << "  line " << ins.line << ": ";
      if (ins.op==qasm::Op::CUSTOM) std::cout << ins.name;
      else std::cout << (int)ins.op << "(" << op_name(ins.op) << ")";
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
  
  static std::string read_all_stdin() {
    std::ostringstream oss;
    oss << std::cin.rdbuf();
    return oss.str();
  }
  
  static std::string read_all_file(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Failed to open: " + path);
    std::ostringstream oss; oss << ifs.rdbuf();
    return oss.str();
  }
  
} // namespace qasm
