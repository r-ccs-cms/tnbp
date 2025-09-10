// main.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "qasm/any.h"
#include "qasm/utility.h"

int main(int argc, char** argv){
  try {
    std::string src;
    if (argc >= 2) {
      src = qasm::read_all_file(argv[1]);     // readin from a file if the filename is written as an argument
    } else {
      src = qasm::read_all_stdin();           // readin from stdin if the filename is not written
    }

    qasm::Program p = qasm::parse_any(src);
    qasm::dump(p);
  } catch(const qasm::ParseError& e){
    std::cerr << e.what() << "\n";
    return 2;
  } catch(const std::exception& e){
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
