// qasm/v2/parser.hpp
#pragma once
#include "../common/reader.hpp"
#include <unordered_map>

namespace qasm::v2 {

class Parser {
 public:
  qasm::Program parse(std::string_view src) {
    r.reset(src);
    r.skipSpacesAndComments();
    if (r.peekWordEq("OPENQASM")) {
      r.readWord();  // OPENQASM
      r.consumeVersionNumberToken();
      r.expect(';');
    }
    prog = {};
    reg_q_sizes.clear(); reg_c_sizes.clear();

    while (!r.eof()) {
      r.skipSpacesAndComments();
      if (r.eof()) break;
      if (r.peekWordEq("include")) { readInclude(); continue; }
      if (r.peekWordEq("qreg"))    { readQReg();    continue; }
      if (r.peekWordEq("creg"))    { readCReg();    continue; }
      if (r.peekWordEq("barrier")) { readBarrier(); continue; }
      if (r.peekWordEq("reset"))   { readReset();   continue; }
      if (r.peekWordEq("measure")) { readMeasureArrow(); continue; }
      readGateOrCustom();
    }
    return prog;
  }

 private:
  // ---- helpers ----
  void readInclude() {
    r.readWord(); // include
    r.skipSpacesAndComments(); r.expect('"');
    while (!r.eof() && r.get()!='"') {}
    r.expect(';');
  }

  void readQReg() {
    r.readWord(); // qreg
    std::string nm = r.readWord();
    r.skipSpacesAndComments(); r.expect('[');
    std::size_t n = (std::size_t)r.parseNumber();
    r.expect(']'); r.expect(';');
    reg_q_sizes[nm]=n; prog.qregs.push_back({nm,n});
  }

  void readCReg() {
    r.readWord(); // creg
    std::string nm = r.readWord();
    r.skipSpacesAndComments(); r.expect('[');
    std::size_t n = (std::size_t)r.parseNumber();
    r.expect(']'); r.expect(';');
    reg_c_sizes[nm]=n; prog.cregs.push_back({nm,n});
  }

  void checkQ(const std::string& reg, std::size_t idx){
    auto it=reg_q_sizes.find(reg); if (it==reg_q_sizes.end()) r.err("unknown qreg: "+reg);
    if (idx>=it->second) r.err("qreg index out of range: "+reg+"["+std::to_string(idx)+"]");
  }
  void checkC(const std::string& reg, std::size_t idx){
    auto it=reg_c_sizes.find(reg); if (it==reg_c_sizes.end()) r.err("unknown creg: "+reg);
    if (idx>=it->second) r.err("creg index out of range: "+reg+"["+std::to_string(idx)+"]");
  }

  std::vector<qasm::QArg> readQArgs(bool allow_many=true){
    std::vector<qasm::QArg> out;
    auto readOne=[&](){
      r.skipSpacesAndComments();
      std::string reg = r.readWord();
      r.skipSpacesAndComments();
      if (r.peek()=='['){
        r.get();
        std::size_t idx = (std::size_t)r.parseNumber();
        r.expect(']');
        checkQ(reg, idx);
        return std::vector<qasm::QArg>{{reg,idx}};
      } else {
        auto it=reg_q_sizes.find(reg); if (it==reg_q_sizes.end()) r.err("unknown qreg: "+reg);
        std::vector<qasm::QArg> v; v.reserve(it->second);
        for (std::size_t k=0;k<it->second;++k) v.push_back({reg,k});
        return v;
      }
    };
    auto first = readOne(); out.insert(out.end(), first.begin(), first.end());
    if (!allow_many) return out;
    r.skipSpacesAndComments();
    while (r.peek()==','){ r.get(); auto nx=readOne(); out.insert(out.end(), nx.begin(), nx.end()); r.skipSpacesAndComments(); }
    return out;
  }

  std::pair<std::string,std::size_t> readCArg(){
    r.skipSpacesAndComments();
    std::string reg = r.readWord(); r.skipSpacesAndComments(); r.expect('[');
    std::size_t idx = (std::size_t)r.parseNumber(); r.expect(']');
    checkC(reg, idx); return {reg, idx};
  }

  void readBarrier() {
    std::size_t L = r.line; r.readWord();
    auto qargs = readQArgs(true); r.expect(';');
    qasm::Instruction in; in.op=qasm::Op::BARRIER; in.qubits=std::move(qargs); in.line=L; in.raw=r.currentRawLine(L);
    prog.instructions.push_back(std::move(in));
  }

  void readReset() {
    std::size_t L = r.line; r.readWord();
    auto qargs = readQArgs(true); r.expect(';');
    for (auto &qa: qargs) {
      qasm::Instruction in; in.op=qasm::Op::RESET; in.qubits={qa}; in.line=L; in.raw=r.currentRawLine(L);
      prog.instructions.push_back(std::move(in));
    }
  }

  void readMeasureArrow() {
    std::size_t L = r.line; r.readWord(); // measure
    auto q = readQArgs(false);
    r.skipSpacesAndComments();
    if (!(r.peek()=='-' && r.i+1<r.src.size() && r.src[r.i+1]=='>')) r.err("expected '->' after measure");
    r.get(); r.get();
    auto [cr,ci]=readCArg(); r.expect(';');
    if (q.size()!=1) r.err("measure expects a single qubit");
    qasm::Instruction in; in.op=qasm::Op::MEASURE; in.qubits={q[0]}; in.has_classical=true; in.c_reg=cr; in.c_index=ci; in.line=L; in.raw=r.currentRawLine(L);
    prog.instructions.push_back(std::move(in));
  }

  qasm::Op mapGate(const std::string& g){
    using R=qasm::common::Reader;
    auto L=[&](const char* s){ return R::eqNoCase(g,s); };
    if (L("u3")) return qasm::Op::U3;
    if (L("u2")) return qasm::Op::U2;
    if (L("u1")) return qasm::Op::U1;
    if (L("rx")) return qasm::Op::RX;
    if (L("ry")) return qasm::Op::RY;
    if (L("rz")) return qasm::Op::RZ;
    if (L("h"))  return qasm::Op::H;
    if (L("x"))  return qasm::Op::X;
    if (L("y"))  return qasm::Op::Y;
    if (L("z"))  return qasm::Op::Z;
    if (L("s"))  return qasm::Op::S;
    if (L("sdg"))return qasm::Op::SDG;
    if (L("t"))  return qasm::Op::T;
    if (L("tdg"))return qasm::Op::TDG;
    if (L("id")) return qasm::Op::ID;
    if (L("cx")||L("cnot")) return qasm::Op::CX;
    if (L("cz")) return qasm::Op::CZ;
    if (L("swap")) return qasm::Op::SWAP;
    if (L("ccx")||L("toffoli")) return qasm::Op::CCX;
    if (L("cswap")) return qasm::Op::CSWAP;
    return qasm::Op::CUSTOM;
  }
  std::size_t gateArity(qasm::Op op){
    switch(op){
      case qasm::Op::U3: case qasm::Op::U2: case qasm::Op::U1:
      case qasm::Op::RX: case qasm::Op::RY: case qasm::Op::RZ:
      case qasm::Op::H: case qasm::Op::X: case qasm::Op::Y: case qasm::Op::Z:
      case qasm::Op::S: case qasm::Op::SDG: case qasm::Op::T: case qasm::Op::TDG:
      case qasm::Op::ID: return 1;
      case qasm::Op::CX: case qasm::Op::CZ: case qasm::Op::SWAP: return 2;
      case qasm::Op::CCX: case qasm::Op::CSWAP: return 3;
      default: return 0;
    }
  }

  void readGateOrCustom(){
    std::size_t L=r.line;
    std::string gate = r.readWord();
    std::vector<double> params;
    r.skipSpacesAndComments();
    if (r.peek()=='(') params = r.parseParamList();
    r.skipSpacesAndComments();
    auto qargs = readQArgs(true);
    r.expect(';');

    qasm::Op op = mapGate(gate);
    std::size_t ar = gateArity(op);

    if (ar>0 && qargs.size()%ar!=0) r.err("gate '"+gate+"' expects multiple of "+std::to_string(ar)+" qubits");

    if (ar==0){
      qasm::Instruction in; in.op=op; if (op==qasm::Op::CUSTOM) in.name=gate; in.params=params; in.qubits=std::move(qargs); in.line=L; in.raw=r.currentRawLine(L);
      prog.instructions.push_back(std::move(in));
    } else {
      for (std::size_t k=0;k<qargs.size();k+=ar){
        qasm::Instruction in; in.op=op; if (op==qasm::Op::CUSTOM) in.name=gate; in.params=params;
        in.qubits.assign(qargs.begin()+k, qargs.begin()+k+ar); in.line=L; in.raw=r.currentRawLine(L);
        prog.instructions.push_back(std::move(in));
      }
    }
  }

 private:
  qasm::common::Reader r;
  qasm::Program prog;
  std::unordered_map<std::string,std::size_t> reg_q_sizes, reg_c_sizes;
};

} // namespace qasm::v2
