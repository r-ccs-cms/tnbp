// qasm/common/reader.hpp
#pragma once
#include "../ir.hpp"
#include <string_view>
#include <string>
#include <vector>
#include <cctype>
#include <cstdlib>
#include <sstream>

namespace qasm::common {

struct Reader {
  std::string_view src{};
  std::size_t i{0}, line{1}, col{1};
  std::vector<std::size_t> raw_line_starts;

  void reset(std::string_view s) {
    src = s; i = 0; line = 1; col = 1;
    raw_line_starts.clear(); raw_line_starts.push_back(0);
  }
  bool eof() const { return i >= src.size(); }
  char peek() const { return eof()? '\0' : src[i]; }
  char get() {
    if (eof()) return '\0';
    char c = src[i++];
    if (c=='\n') { line++; col=1; raw_line_starts.push_back(i); }
    else col++;
    return c;
  }
  void unget() { if (i>0) { --i; if (col>0) --col; } }

  static bool isIdentStart(char c){ return std::isalpha((unsigned char)c) || c=='_'; }
  static bool isIdentChar (char c){ return std::isalnum((unsigned char)c) || c=='_' || c=='$'; }

  void skipSpaces() { while (!eof() && std::isspace((unsigned char)peek())) get(); }
  void skipSpacesAndComments() {
    while (true) {
      skipSpaces();
      if (peek()=='/' && i+1<src.size() && src[i+1]=='/') { while (!eof() && get()!='\n'){} continue; }
      if (peek()=='/' && i+1<src.size() && src[i+1]=='*') {
        get(); get();
        while (!eof()) { char c=get(); if (c=='*' && peek()=='/') { get(); break; } }
        continue;
      }
      break;
    }
  }

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
    std::size_t si=i, sl=line, sc=col;
    std::string t = readWord();
    bool ok = eqNoCase(t, w);
    i=si; line=sl; col=sc;
    return ok;
  }

  void expect(char ch) {
    skipSpacesAndComments();
    if (peek()!=ch) { std::string m="expected '"; m.push_back(ch); m+="'"; err(m); }
    get();
  }

  // Very permissive OpenQASM header version eater (e.g. "2.0" or "3")
  void consumeVersionNumberToken() {
    skipSpacesAndComments();
    if (peek()=='"') { // not expected but be safe
      get(); while(!eof() && get()!='"'){}
      return;
    }
    while (!eof()) {
      char c=peek();
      if (std::isdigit((unsigned char)c) || c=='.') get();
      else break;
    }
  }

  double parseExpr() { double v=parseTerm(); while(true){ skipSpacesAndComments(); char c=peek();
    if (c=='+'){ get(); v+=parseTerm(); } else if (c=='-'){ get(); v-=parseTerm(); } else break; } return v; }
  double parseTerm() { double v=parseFactor(); while(true){ skipSpacesAndComments(); char c=peek();
    if (c=='*'){ get(); v*=parseFactor(); } else if (c=='/'){ get(); double d=parseFactor(); if (d==0.0) err("division by zero"); v/=d; } else break; } return v; }

  double parseNumber() {
    skipSpacesAndComments();
    std::size_t s=i; bool has_dot=false, has_exp=false;
    if (peek()=='+'||peek()=='-') get();
    while(!eof()){
      char c=peek();
      if (std::isdigit((unsigned char)c)) get();
      else if (c=='.' && !has_dot) { has_dot=true; get(); }
      else if ((c=='e'||c=='E') && !has_exp) { has_exp=true; get(); if (peek()=='+'||peek()=='-') get(); }
      else break;
    }
    if (i==s) err("number expected");
    return std::strtod(std::string(src.substr(s,i-s)).c_str(), nullptr);
  }

  double parseFactor() {
    skipSpacesAndComments();
    char c=peek();
    if (c=='('){ get(); double v=parseExpr(); expect(')'); return v; }
    if (c=='+'||c=='-'){ bool neg=(c=='-'); get(); double v=parseFactor(); return neg? -v : v; }
    if (isIdentStart(c)){
      std::string w = readWord();
      if (eqNoCase(w,"pi"))  return 3.14159265358979323846;
      if (eqNoCase(w,"tau")) return 6.28318530717958647692;
      err("unknown symbol in expression: "+w);
    }
    return parseNumber();
  }

  std::vector<double> parseParamList() {
    skipSpacesAndComments();
    std::vector<double> ps; expect('('); skipSpacesAndComments();
    if (peek()!=')'){ while(true){ ps.push_back(parseExpr()); skipSpacesAndComments(); if (peek()==','){ get(); continue; } break; } }
    expect(')'); return ps;
  }

  std::string currentRawLine(std::size_t line_no) const {
    if (raw_line_starts.empty() || line_no==0 || line_no>raw_line_starts.size()) return {};
    std::size_t start = raw_line_starts[line_no-1];
    std::size_t end = (line_no<raw_line_starts.size()? raw_line_starts[line_no]-1 : src.size());
    while (end<src.size() && src[end]!='\n') end++;
    return std::string(src.substr(start, end-start));
  }

  [[noreturn]] void err(const std::string& msg) const {
    throw qasm::ParseError(line, col, msg);
  }

  static bool eqNoCase(const std::string& a, std::string_view b){
    if (a.size()!=b.size()) return false;
    for (std::size_t k=0;k<a.size();++k)
      if (std::tolower((unsigned char)a[k]) != std::tolower((unsigned char)b[k])) return false;
    return true;
  }
};

} // namespace qasm::common
