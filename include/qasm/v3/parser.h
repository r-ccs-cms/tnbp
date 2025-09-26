// qasm/v3/parser.hpp
#pragma once
#include "../common/reader.h"
#include <unordered_map>

namespace qasm::v3 {

class Parser {
 public:
  qasm::Program parse(std::string_view src) {
    r.reset(src);
    r.skipSpacesAndComments();
    if (!r.peekWordEq("OPENQASM")) r.err("OPENQASM header required for v3");
    r.readWord();                   // OPENQASM
    r.consumeVersionNumberToken();  // "3" など
    r.expect(';');

    prog = {};
    reg_q_sizes.clear();
    reg_c_sizes.clear();

    while (!r.eof()) {
      r.skipSpacesAndComments();
      if (r.eof()) break;

      // --- 宣言や単純文 ---
      if (r.peekWordEq("include")) { readInclude(); continue; }
      if (r.peekWordEq("qubit"))   { readQubitDecl(); continue; }
      if (r.peekWordEq("barrier")) { readBarrier();   continue; }
      if (r.peekWordEq("reset"))   { readReset();     continue; }
      if (r.peekWordEq("measure")) { readMeasureArrowCompat(); continue; }

      // ★ 代入測定（bit/creg の LHS）を bit/creg 宣言よりも先に判定
      if (looksLikeBitMeasureAssign() || looksLikeCregMeasureAssign()) {
        readMeasureAssign();
        continue;
      }

      // bit 宣言（通常の宣言）
      if (r.peekWordEq("bit")) { readBitDecl(); continue; }

      // ゲート修飾子は非対応（限定展開を後で入れる場合はここを差し替え）
      {
        std::size_t si=r.i, sl=r.line, sc=r.col;
        std::string w = r.readWord();
        if (ieq(w,"ctrl") || ieq(w,"inv") || ieq(w,"pow")) {
          r.err("gate modifiers (ctrl/inv/pow) are not supported in this subset");
        }
        // 先読みを戻す
        r.i=si; r.line=sl; r.col=sc;
      }

      // ゲート呼び出し
      readGateOrCustom();
    }

    return prog;
  }

 private:
  // ===== ヘルパ：大文字小文字無視比較 =====
  static bool ieq(const std::string& a, std::string_view b){
    return qasm::common::Reader::eqNoCase(a,b);
  }

  // ===== include =====
  void readInclude(){
    r.readWord(); // include
    r.skipSpacesAndComments(); r.expect('"');
    while(!r.eof() && r.get()!='"'){}
    r.expect(';');
  }

  // ===== 宣言（qubit/bit）=====
  void readQubitDecl(){
    r.readWord(); // qubit
    r.skipSpacesAndComments();
    std::size_t size = 1;
    if (r.peek()=='['){ r.get(); size = (std::size_t)r.parseNumber(); r.expect(']'); }
    std::string name = r.readWord();
    r.expect(';');
    reg_q_sizes[name]=size;
    prog.qregs.push_back({name,size});
  }

  void readBitDecl(){
    r.readWord(); // bit
    r.skipSpacesAndComments();
    std::size_t size = 1;
    if (r.peek()=='['){ r.get(); size = (std::size_t)r.parseNumber(); r.expect(']'); }
    std::string name = r.readWord();
    r.expect(';');
    reg_c_sizes[name]=size;
    prog.cregs.push_back({name,size});
  }

  // ===== barrier/reset =====
  void readBarrier(){
    std::size_t L=r.line; r.readWord();
    auto qargs = readQArgs(true);
    r.expect(';');
    qasm::Instruction in;
    in.op = qasm::Op::BARRIER;
    in.qubits = std::move(qargs);
    in.line = L; in.raw = r.currentRawLine(L);
    prog.instructions.push_back(std::move(in));
  }

  void readReset(){
    std::size_t L=r.line; r.readWord();
    auto qargs = readQArgs(true);
    r.expect(';');
    for (auto &qa: qargs) {
      qasm::Instruction in;
      in.op = qasm::Op::RESET;
      in.qubits = {qa};
      in.line = L; in.raw = r.currentRawLine(L);
      prog.instructions.push_back(std::move(in));
    }
  }

  // ===== measure 互換 (矢印形) =====
  void readMeasureArrowCompat(){
    std::size_t L=r.line; r.readWord(); // measure
    auto q = readQArgs(false);
    r.skipSpacesAndComments();
    if (!(r.peek()=='-' && r.i+1<r.src.size() && r.src[r.i+1]=='>'))
      r.err("expected '->' after measure");
    r.get(); r.get(); // ->
    auto [cr,ci] = readCArg();
    r.expect(';');
    if (q.size()!=1) r.err("measure expects a single qubit");
    qasm::Instruction in;
    in.op = qasm::Op::MEASURE;
    in.qubits = {q[0]};
    in.has_classical = true;
    in.c_reg = cr; in.c_index = ci;
    in.line = L; in.raw = r.currentRawLine(L);
    prog.instructions.push_back(std::move(in));
  }

  // ===== 代入測定の“見た目”判定 =====
  bool looksLikeBitMeasureAssign(){
    // 例: "bit x = measure q[2];"
    std::size_t si=r.i, sl=r.line, sc=r.col;
    r.skipSpacesAndComments();
    if (!r.peekWordEq("bit")) { r.i=si; r.line=sl; r.col=sc; return false; }
    r.readWord(); // bit
    r.skipSpacesAndComments();
    if (!qasm::common::Reader::isIdentStart(r.peek())) { r.i=si; r.line=sl; r.col=sc; return false; }
    (void)r.readWord(); // name
    r.skipSpacesAndComments();
    bool ok = (r.peek()=='=');
    r.i=si; r.line=sl; r.col=sc;
    return ok;
  }

  bool looksLikeCregMeasureAssign(){
    // 例: "c[i] = measure q[j];" / "c = measure q[0];"（c が size 1 のとき）
    std::size_t si=r.i, sl=r.line, sc=r.col;
    r.skipSpacesAndComments();
    if (!qasm::common::Reader::isIdentStart(r.peek())) { r.i=si; r.line=sl; r.col=sc; return false; }
    (void)r.readWord(); // creg name (仮)
    r.skipSpacesAndComments();
    if (r.peek()=='['){ // index 部
      r.get();
      (void)r.parseNumber();
      r.expect(']');
    }
    r.skipSpacesAndComments();
    bool ok = (r.peek()=='=');
    r.i=si; r.line=sl; r.col=sc;
    return ok;
  }

  // ===== 代入測定の本体 =====
  void readMeasureAssign(){
    // c[i] = measure q[j];  |  bit x = measure q[k];
    std::size_t L=r.line;
    auto [cr, ci] = readCArgOrBitDeclLHS();
    r.skipSpacesAndComments(); r.expect('=');
    r.skipSpacesAndComments();
    if (!r.peekWordEq("measure")) r.err("expected 'measure' in assignment");
    r.readWord(); // measure
    auto q = readQArgs(false);
    r.expect(';');
    if (q.size()!=1) r.err("measure expects a single qubit");
    qasm::Instruction in;
    in.op = qasm::Op::MEASURE;
    in.qubits = {q[0]};
    in.has_classical = true;
    in.c_reg = cr; in.c_index = ci;
    in.line = L; in.raw = r.currentRawLine(L);
    prog.instructions.push_back(std::move(in));
  }

  // LHS: "bit x"（宣言兼ねる）or 既存 "c[i]" / "c" (size==1)
  std::pair<std::string,std::size_t> readCArgOrBitDeclLHS(){
    r.skipSpacesAndComments();
    if (r.peekWordEq("bit")){
      r.readWord(); // bit
      r.skipSpacesAndComments();
      if (r.peek()=='[') r.err("bit array not allowed on LHS with declaration here");
      std::string name = r.readWord();
      if (!reg_c_sizes.count(name)) {
        reg_c_sizes[name]=1;
        prog.cregs.push_back({name,1});
      }
      return {name,0};
    }
    return readCArg(); // 既存 creg
  }

  // c[i] or c (size==1 の省略形)
  std::pair<std::string,std::size_t> readCArg(){
    r.skipSpacesAndComments();
    std::string reg = r.readWord();
    r.skipSpacesAndComments();
    if (r.peek()!='['){
      auto it=reg_c_sizes.find(reg);
      if (it==reg_c_sizes.end()) r.err("unknown bit/creg: "+reg);
      if (it->second!=1) r.err("missing index for creg: "+reg);
      return {reg,0};
    }
    r.get();
    std::size_t idx=(std::size_t)r.parseNumber();
    r.expect(']');
    checkC(reg, idx);
    return {reg, idx};
  }

  // ===== QArgs（単体/複数, スライス [lo:hi) 対応）=====
  std::vector<qasm::QArg> readQArgs(bool allow_many){
    std::vector<qasm::QArg> out;

    auto readOne=[&](){
      r.skipSpacesAndComments();
      std::string reg = r.readWord();
      r.skipSpacesAndComments();
      if (r.peek()=='['){
        r.get();
        // 単一 index か slice [lo:hi)
        r.skipSpacesAndComments();
        std::size_t lo = (std::size_t)r.parseNumber();
        r.skipSpacesAndComments();
        if (r.peek()==':'){
          r.get();
          r.skipSpacesAndComments();
          std::size_t hi = (std::size_t)r.parseNumber();
          r.expect(']');
          auto it=reg_q_sizes.find(reg); if (it==reg_q_sizes.end()) r.err("unknown qreg: "+reg);
          if (hi<lo || hi>it->second) r.err("slice out of range: "+reg+"["+std::to_string(lo)+":"+std::to_string(hi)+"]");
          std::vector<qasm::QArg> v; v.reserve(hi-lo);
          for (std::size_t k=lo;k<hi;++k) v.push_back({reg,k});
          return v;
        } else {
          std::size_t idx = lo;
          r.expect(']');
          checkQ(reg, idx);
          return std::vector<qasm::QArg>{{reg,idx}};
        }
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
    while (r.peek()==','){
      r.get();
      auto nx = readOne();
      out.insert(out.end(), nx.begin(), nx.end());
      r.skipSpacesAndComments();
    }
    return out;
  }

  // ===== gate 呼び出し =====
  void readGateOrCustom(){
    std::size_t L=r.line;
    std::string gate = r.readWord();
    if (ieq(gate,"ctrl") || ieq(gate,"inv") || ieq(gate,"pow")) {
      r.err("gate modifiers are not supported in this subset");
    }
    std::vector<double> params;
    r.skipSpacesAndComments();
    if (r.peek()=='(') params = r.parseParamList();
    r.skipSpacesAndComments();
    auto qargs = readQArgs(true);
    r.expect(';');

    qasm::Op op = mapGate(gate);
    std::size_t ar = gateArity(op);
    if (ar>0 && qargs.size()%ar!=0)
      r.err("gate '"+gate+"' expects multiple of "+std::to_string(ar)+" qubits");

    if (ar==0){
      qasm::Instruction in;
      in.op = op;
      if (op==qasm::Op::CUSTOM) in.name=gate;
      in.params = params;
      in.qubits = std::move(qargs);
      in.line = L; in.raw = r.currentRawLine(L);
      prog.instructions.push_back(std::move(in));
    } else {
      for (std::size_t k=0;k<qargs.size();k+=ar){
        qasm::Instruction in;
        in.op = op;
        if (op==qasm::Op::CUSTOM) in.name=gate;
        in.params = params;
        in.qubits.assign(qargs.begin()+k, qargs.begin()+k+ar);
        in.line = L; in.raw = r.currentRawLine(L);
        prog.instructions.push_back(std::move(in));
      }
    }
  }

  // ===== gate マッピング =====
  qasm::Op mapGate(const std::string& g){
    using R=qasm::common::Reader;
    auto L=[&](const char* s){ return R::eqNoCase(g,s); };
    if (L("u3") || L("u")) return qasm::Op::U3;     // In Qiskit u = U3
    if (L("u2"))           return qasm::Op::U2;
    if (L("u1") || L("p")) return qasm::Op::U1;   // In Qiskit p(t) = u1(t)
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

  // ===== index 範囲チェック =====
  void checkQ(const std::string& reg, std::size_t idx){
    auto it=reg_q_sizes.find(reg);
    if (it==reg_q_sizes.end()) r.err("unknown qreg: "+reg);
    if (idx>=it->second) r.err("qreg index out of range: "+reg+"["+std::to_string(idx)+"]");
  }

  void checkC(const std::string& reg, std::size_t idx){
    auto it=reg_c_sizes.find(reg);
    if (it==reg_c_sizes.end()) r.err("unknown creg/bit: "+reg);
    if (idx>=it->second) r.err("creg index out of range: "+reg+"["+std::to_string(idx)+"]");
  }

 private:
  qasm::common::Reader r;
  qasm::Program prog;
  std::unordered_map<std::string,std::size_t> reg_q_sizes, reg_c_sizes;
};

} // namespace qasm::v3
