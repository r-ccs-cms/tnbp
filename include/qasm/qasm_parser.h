// ============================================================
// Coverage & Design Notes
// ============================================================
//
// * Expansion Policy:
//   Statements such as `h q;` that apply to an entire register
//   are expanded into individual instructions for each qubit,
//   making the resulting IR easier to process.
//
// * Angle Expressions:
//   Supports arithmetic (+, -, *, /), parentheses, unary ±,
//   and symbolic constants `pi` and `tau`. Example:
//   `-(pi/2) + 0.1`. Any unknown symbol raises an error.
//
// * Custom Gates:
//   Unknown gate names are treated as `Op::CUSTOM`.
//   Their arity is set to zero (free), leaving interpretation
//   up to the user. For stricter handling, extend `gateArity()`
//   or add gate-definition parsing.
//
// * Error Messages:
//   Errors throw `ParseError(line:col, message)`
//   so the failing position in the QASM source can be located.
//
// * Register-Wide Application:
//   When a register is specified without an index (e.g. `q`),
//   the parser expands it into per-qubit arguments automatically.
//
// * Dependencies:
//   Header-only, using only the C++17 standard library.
//   No external libraries required.
//
// * Debugging Aid:
//   Each instruction keeps its original line number and raw
//   source line for easier debugging.
//
// * Extension Ideas:
//   - Support for `gate` definitions (`gate name(params) qargs { ... }`)
//   - Parameter identifiers (variables like `theta`)
//   - Conversion utilities from `Program` to backend IR
//   - Experimental support for OpenQASM 3 features
//
// ============================================================
#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <optional>
#include <limits>

namespace qasm {

struct QReg { std::string name; std::size_t size{}; };
struct CReg { std::string name; std::size_t size{}; };

enum class Op {
  U3,U2,U1, RX,RY,RZ, H,X,Y,Z,S,SDG,T,TDG, CX,CZ, SWAP, CCX, CSWAP, ID,
  MEASURE, BARRIER, RESET, CUSTOM
};

struct QArg {
  std::string reg;
  std::size_t index{};
};

struct Instruction {
  Op op{};
  std::string name;                 // gate name when CUSTOM 
  std::vector<double> params;       // parameters
  std::vector<QArg> qubits;         // qubits where gate acts
  bool has_classical{false};        // when it is MEASURE
  std::string c_reg;                // when it is MEASURE
  std::size_t c_index{};
  std::size_t line{};               // line number of original source
  std::string raw;                  // original line for debug
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

class Parser {
 public:
  Program parse(std::string_view src) {
    reset(src);
    skipSpacesAndComments();
    if (peekWordEq("OPENQASM")) {
      readWord();  // OPENQASM
      expectNumberOrVersion();
      expect(';');
    }
    // skip include (qelib1.inc etc.)
    while (!eof()) {
      skipSpacesAndComments();
      if (eof()) break;
      if (peekWordEq("include")) { readInclude(); continue; }
      if (peekWordEq("qreg"))    { readQReg();    continue; }
      if (peekWordEq("creg"))    { readCReg();    continue; }
      if (peekWordEq("barrier")) { readBarrier(); continue; }
      if (peekWordEq("reset"))   { readReset();   continue; }
      if (peekWordEq("measure")) { readMeasure(); continue; }
      // ゲート適用
      readGateOrCustom();
    }
    return prog_;
  }

  Program parseFile(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Failed to open: "+path);
    std::ostringstream oss; oss << ifs.rdbuf();
    return parse(oss.str());
  }

 private:
  // ===== lexer-like helpers =====
  void reset(std::string_view src) {
    src_ = src;
    i_ = 0; line_ = 1; col_ = 1;
    prog_ = Program{};
    raw_line_starts_.clear();
    raw_line_starts_.push_back(0);
  }
  bool eof() const { return i_ >= src_.size(); }
  char peek() const { return eof()? '\0' : src_[i_]; }
  char get() {
    if (eof()) return '\0';
    char c = src_[i_++];
    if (c=='\n') { line_++; col_=1; raw_line_starts_.push_back(i_); }
    else col_++;
    return c;
  }
  void unget() { // only for 1-char step back
    if (i_==0) return;
    i_--;
    if (i_>0 && src_[i_]=='\n') { /* shouldn't happen for our usage */ }
    // simplistic col_ tracking is fine (we only use for errors)
  }

  void skipSpaces() { while (!eof() && std::isspace((unsigned char)peek())) get(); }

  void skipSpacesAndComments() {
    while (true) {
      skipSpaces();
      if (peek()=='/' && i_+1<src_.size() && src_[i_+1]=='/') {
        // line comment
        while (!eof() && get()!='\n') {}
        continue;
      }
      if (peek()=='/' && i_+1<src_.size() && src_[i_+1]=='*') {
        // block comment
        get(); get();
        while (!eof()) {
          char c = get();
          if (c=='*' && peek()=='/') { get(); break; }
        }
        continue;
      }
      break;
    }
  }

  bool isIdentStart(char c) const { return std::isalpha((unsigned char)c) || c=='_'; }
  bool isIdentChar(char c)  const { return std::isalnum((unsigned char)c) || c=='_' || c=='$'; }

  std::string readWord() {
    skipSpacesAndComments();
    if (!isIdentStart(peek())) err("identifier expected");
    std::string w; w.push_back(get());
    while (isIdentChar(peek())) w.push_back(get());
    return w;
  }

  bool peekWordEq(std::string_view w) {
    skipSpacesAndComments();
    if (!isIdentStart(peek())) return false;
    std::size_t save_i = i_, save_line = line_, save_col = col_;
    std::string t = readWord();
    bool ok = eqNoCase(t, w);
    i_ = save_i; line_ = save_line; col_ = save_col;
    return ok;
  }

  bool eqNoCase(const std::string& a, std::string_view b) const {
    if (a.size()!=b.size()) return false;
    for (std::size_t k=0;k<a.size();++k) {
      if (std::tolower((unsigned char)a[k]) != std::tolower((unsigned char)b[k])) return false;
    }
    return true;
  }

  void expect(char ch) {
    skipSpacesAndComments();
    if (peek()!=ch) {
      std::string m = "expected '"; m.push_back(ch); m += "'";
      err(m);
    }
    get();
  }

  void expectNumberOrVersion() {
    // consumes something like "2.0" or a number token
    skipSpacesAndComments();
    // very permissive: sequence [0-9|.|a-z] until ';' or space
    while (!eof()) {
      char c = peek();
      if (std::isdigit((unsigned char)c) || c=='.') { get(); }
      else break;
    }
  }

  // ===== expression parser (for angles) =====
  double parseExpr() { // +,-
    double v = parseTerm();
    while (true) {
      skipSpacesAndComments();
      char c = peek();
      if (c=='+') { get(); v += parseTerm(); }
      else if (c=='-') { get(); v -= parseTerm(); }
      else break;
    }
    return v;
  }
  double parseTerm() { // *,/
    double v = parseFactor();
    while (true) {
      skipSpacesAndComments();
      char c = peek();
      if (c=='*') { get(); v *= parseFactor(); }
      else if (c=='/') { get(); v /= parseFactorSafe(); }
      else break;
    }
    return v;
  }
  double parseFactorSafe() {
    double d = parseFactor();
    if (d==0.0) err("division by zero in angle expression");
    return d;
  }
  double parseNumber() {
    skipSpacesAndComments();
    std::size_t start = i_;
    bool has_dot=false, has_exp=false;
    if (peek()=='+' || peek()=='-') get();
    while (!eof()) {
      char c = peek();
      if (std::isdigit((unsigned char)c)) { get(); }
      else if (c=='.' && !has_dot) { has_dot=true; get(); }
      else if ((c=='e' || c=='E') && !has_exp) {
        has_exp=true; get();
        if (peek()=='+' || peek()=='-') get();
      } else break;
    }
    if (i_==start) err("number expected");
    double v = std::strtod(std::string(src_.substr(start, i_-start)).c_str(), nullptr);
    return v;
  }
  double parseFactor() {
    skipSpacesAndComments();
    char c = peek();
    if (c=='(') {
      get();
      double v = parseExpr();
      expect(')');
      return v;
    }
    if (c=='+' || c=='-') { // unary
      bool neg = (c=='-'); get();
      double v = parseFactor();
      return neg? -v : v;
    }
    if (isIdentStart(c)) {
      std::string w = readWord();
      if (eqNoCase(w, "pi"))  return 3.141592653589793238462643383279502884;
      if (eqNoCase(w, "tau")) return 6.283185307179586476925286766559005768; // 2*pi
      // not implemented for unkown parameters -> error
      err("unknown symbol in expression: "+w);
    }
    return parseNumber();
  }

  std::vector<double> parseParamList() {
    skipSpacesAndComments();
    std::vector<double> ps;
    expect('(');
    skipSpacesAndComments();
    if (peek()!=')') {
      while (true) {
        double v = parseExpr();
        ps.push_back(v);
        skipSpacesAndComments();
        if (peek()==',') { get(); continue; }
        break;
      }
    }
    expect(')');
    return ps;
  }

  // ===== QASM constructs =====
  void readInclude() {
    readWord(); // include
    skipSpacesAndComments();
    expect('"');
    // skip until next quote
    while (!eof() && get()!='"') {}
    expect(';');
  }

  void readQReg() {
    readWord(); // qreg
    std::string nm = readWord();
    skipSpacesAndComments();
    expect('[');
    std::size_t n = (std::size_t)parseNumber(); // only integer expected
    expect(']');
    expect(';');
    reg_q_sizes_[nm] = n;
    prog_.qregs.push_back(QReg{nm,n});
  }

  void readCReg() {
    readWord(); // creg
    std::string nm = readWord();
    skipSpacesAndComments();
    expect('[');
    std::size_t n = (std::size_t)parseNumber();
    expect(']');
    expect(';');
    reg_c_sizes_[nm] = n;
    prog_.cregs.push_back(CReg{nm,n});
  }

  std::vector<QArg> readQArgListOneOrMany(bool allow_many=true) {
    // q[0] | q  | q[1],r[2]...
    std::vector<QArg> out;
    auto readOne = [&]() -> std::vector<QArg> {
      skipSpacesAndComments();
      std::string reg = readWord();
      skipSpacesAndComments();
      if (peek()=='[') {
        get();
        std::size_t idx = (std::size_t)parseNumber();
        expect(']');
        checkQIndex(reg, idx);
        return { QArg{reg, idx} };
      } else {
        // whole register
        auto it = reg_q_sizes_.find(reg);
        if (it==reg_q_sizes_.end()) err("unknown qreg: "+reg);
        std::vector<QArg> many; many.reserve(it->second);
        for (std::size_t k=0;k<it->second;++k) many.push_back(QArg{reg,k});
        return many;
      }
    };

    auto first = readOne();
    out.insert(out.end(), first.begin(), first.end());
    if (!allow_many) return out;

    skipSpacesAndComments();
    while (peek()==',') {
      get();
      auto next = readOne();
      out.insert(out.end(), next.begin(), next.end());
      skipSpacesAndComments();
    }
    return out;
  }

  std::pair<std::string,std::size_t> readCArg() {
    // c[0]
    skipSpacesAndComments();
    std::string reg = readWord();
    skipSpacesAndComments();
    expect('[');
    std::size_t idx = (std::size_t)parseNumber();
    expect(']');
    checkCIndex(reg, idx);
    return {reg, idx};
  }

  void readBarrier() {
    std::size_t start_line = line_;
    readWord(); // barrier
    auto qargs = readQArgListOneOrMany(/*allow_many=*/true);
    expect(';');
    Instruction in;
    in.op = Op::BARRIER;
    in.qubits = std::move(qargs);
    in.line = start_line;
    in.raw  = currentRawLine(start_line);
    prog_.instructions.push_back(std::move(in));
  }

  void readReset() {
    std::size_t start_line = line_;
    readWord(); // reset
    auto qargs = readQArgListOneOrMany(/*allow_many=*/true);
    expect(';');
    for (auto &qa : qargs) {
      Instruction in;
      in.op = Op::RESET;
      in.qubits = {qa};
      in.line = start_line;
      in.raw  = currentRawLine(start_line);
      prog_.instructions.push_back(std::move(in));
    }
  }

  void readMeasure() {
    std::size_t start_line = line_;
    readWord(); // measure
    auto q = readQArgListOneOrMany(/*allow_many=*/false); // it is only one in the case of OpenQASM 2.0
    skipSpacesAndComments();
    if (!(peek()=='-' && i_+1<src_.size() && src_[i_+1]=='>')) err("expected '->' after measure");
    get(); get(); // ->
    auto [cr, ci] = readCArg();
    expect(';');
    if (q.size()!=1) err("measure expects a single qubit (use q[i])");
    Instruction in;
    in.op = Op::MEASURE;
    in.qubits = {q[0]};
    in.has_classical = true;
    in.c_reg = cr;
    in.c_index = ci;
    in.line = start_line;
    in.raw  = currentRawLine(start_line);
    prog_.instructions.push_back(std::move(in));
  }

  void readGateOrCustom() {
    std::size_t start_line = line_;
    std::string gate = readWord(); // name of gate
    std::vector<double> params;
    skipSpacesAndComments();
    if (peek()=='(') params = parseParamList();
    skipSpacesAndComments();
    auto qargs = readQArgListOneOrMany(/*allow_many=*/true);
    expect(';');

    // it is already implemented for each registart (see inside of readQArgListOneOrMany)

    // check compatibility of gates which has more than 2 qargs and integrate Instruction
    Op op = mapGate(gate);
    std::size_t arity = gateArity(op, gate);

    if (arity>0 && qargs.size()%arity!=0) {
      err("gate '"+gate+"' expects qubit count multiple of "+std::to_string(arity));
    }

    // introduce separation in accordance to arity
    if (arity==0) { // arbitrary number is allows (e.g. barrier is already taken into account)
      Instruction in;
      in.op = op;
      in.name = (op==Op::CUSTOM? gate : "");
      in.params = std::move(params);
      in.qubits = std::move(qargs);
      in.line = start_line;
      in.raw  = currentRawLine(start_line);
      prog_.instructions.push_back(std::move(in));
    } else {
      for (std::size_t k=0;k<qargs.size();k+=arity) {
        Instruction in;
        in.op = op;
        in.name = (op==Op::CUSTOM? gate : "");
        in.params = params;
        in.qubits.assign(qargs.begin()+k, qargs.begin()+k+arity);
        in.line = start_line;
        in.raw  = currentRawLine(start_line);
        prog_.instructions.push_back(std::move(in));
      }
    }
  }

  // ===== utilities =====
  Op mapGate(const std::string& g) {
    auto L = [&](const char* s){ return eqNoCase(g,s); };
    if (L("u3")) return Op::U3;
    if (L("u2")) return Op::U2;
    if (L("u1")) return Op::U1;
    if (L("rx")) return Op::RX;
    if (L("ry")) return Op::RY;
    if (L("rz")) return Op::RZ;
    if (L("h"))  return Op::H;
    if (L("x"))  return Op::X;
    if (L("y"))  return Op::Y;
    if (L("z"))  return Op::Z;
    if (L("s"))  return Op::S;
    if (L("sdg"))return Op::SDG;
    if (L("t"))  return Op::T;
    if (L("tdg"))return Op::TDG;
    if (L("id")) return Op::ID;
    if (L("cx")||L("cnot")) return Op::CX;
    if (L("cz")) return Op::CZ;
    if (L("swap")) return Op::SWAP;
    if (L("ccx")||L("toffoli")) return Op::CCX;
    if (L("cswap")) return Op::CSWAP;
    // unknown is CUSTOM (number of parameters and arity is interpreted by user)
    return Op::CUSTOM;
  }

  std::size_t gateArity(Op op, const std::string& name) {
    switch (op) {
      case Op::U3: case Op::U2: case Op::U1:
      case Op::RX: case Op::RY: case Op::RZ:
      case Op::H: case Op::X: case Op::Y: case Op::Z:
      case Op::S: case Op::SDG: case Op::T: case Op::TDG:
      case Op::ID:
        return 1;
      case Op::CX: case Op::CZ: case Op::SWAP:
        return 2;
      case Op::CCX: case Op::CSWAP:
        return 3;
      case Op::MEASURE: case Op::BARRIER: case Op::RESET:
        return 0; // We assume it is not come to here
      case Op::CUSTOM:
        return 0;
    }
    return 0;
  }

  void checkQIndex(const std::string& r, std::size_t idx) {
    auto it = reg_q_sizes_.find(r);
    if (it==reg_q_sizes_.end()) err("unknown qreg: "+r);
    if (idx >= it->second) err("qreg index out of range: "+r+"["+std::to_string(idx)+"]");
  }
  void checkCIndex(const std::string& r, std::size_t idx) {
    auto it = reg_c_sizes_.find(r);
    if (it==reg_c_sizes_.end()) err("unknown creg: "+r);
    if (idx >= it->second) err("creg index out of range: "+r+"["+std::to_string(idx)+"]");
  }

  [[noreturn]] void err(const std::string& msg) const {
    throw ParseError(line_, col_, msg);
  }

  std::string currentRawLine(std::size_t line_no) const {
    if (raw_line_starts_.empty() || line_no==0 || line_no>raw_line_starts_.size())
      return {};
    std::size_t start = raw_line_starts_[line_no-1];
    std::size_t end = (line_no<raw_line_starts_.size()? raw_line_starts_[line_no]-1 : src_.size());
    // extend to next newline purely for debug readability
    while (end<src_.size() && src_[end]!='\n') end++;
    return std::string(src_.substr(start, end-start));
  }

  // ===== state =====
  std::string_view src_;
  std::size_t i_{0};
  std::size_t line_{1}, col_{1};
  Program prog_;
  std::unordered_map<std::string,std::size_t> reg_q_sizes_;
  std::unordered_map<std::string,std::size_t> reg_c_sizes_;
  std::vector<std::size_t> raw_line_starts_;
};

} // namespace qasm
