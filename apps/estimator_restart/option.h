#pragma once
#include <charconv>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Option {

  // inputs for problem settings
  std::string ownerfile;
  std::string backend = "default";
  std::string circuit = "circuit.qasm";
  std::string sparsepauli = "sparsepauliop.txt";

  // io and variables for restart
  std::string portable_save, portable_load;
  bool step_start_set = false;
  std::string savename;
  std::string loadname;
  size_t step_start = 0;
  size_t step_end = 0;

  // inputs for controlled parameters
  std::vector<int> measurement_barrier;
  size_t max_bp_iterations = 50;
  size_t max_bond_dim = 100;
  double bp_tolerance = 1.0e-8;
  double sv_min = 1.0e-8;
  double truncation_error = 1.0e-8;
};

Option generate_options(int argc, char *argv[]) {

  Option option;
  auto count = [](const std::string &flag, const std::string &value) {
    std::size_t n = 0;
    auto parsed = std::from_chars(value.data(), value.data() + value.size(), n);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size())
      throw std::invalid_argument(flag + " requires a nonnegative integer");
    return n;
  };
  auto real = [](const std::string &flag, const std::string &value) {
    std::size_t used = 0;
    auto n = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(n) || n < 0)
      throw std::invalid_argument(flag +
                                  " requires a finite nonnegative number");
    return n;
  };
  for (int i = 1; i < argc; ++i) {
    const std::string flag = argv[i];
    if (i + 1 >= argc || std::string(argv[i + 1]).rfind("--", 0) == 0)
      throw std::invalid_argument(flag + " requires a value");
    const std::string value = argv[++i];
    if (value.empty())
      throw std::invalid_argument(flag + " requires a nonempty value");
    if (flag == "--ownerfile")
      option.ownerfile=value;
    else if (flag == "--backend")
      option.backend = value;
    else if (flag == "--circuit")
      option.circuit = value;
    else if (flag == "--sparse_pauli")
      option.sparsepauli = value;
    else if (flag == "--savename")
      option.savename = value;
    else if (flag == "--loadname")
      option.loadname = value;
    else if (flag == "--portable_save")
      option.portable_save = value;
    else if (flag == "--portable_load")
      option.portable_load = value;
    else if (flag == "--step_start") {
      option.step_start_set = true;
      option.step_start = count(flag, value);
    } else if (flag == "--step_end")
      option.step_end = count(flag, value);
    else if (flag == "--max_bp_iterations")
      option.max_bp_iterations = count(flag, value);
    else if (flag == "--max_bond_dim") {
      option.max_bond_dim = count(flag, value);
      if (option.max_bond_dim == 0 ||
          option.max_bond_dim > std::numeric_limits<int>::max())
        throw std::invalid_argument(flag + " exceeds bond dimension range");
    } else if (flag == "--bp_tolerance")
      option.bp_tolerance = real(flag, value);
    else if (flag == "--sv_min")
      option.sv_min = real(flag, value);
    else if (flag == "--truncation_error")
      option.truncation_error = real(flag, value);
    else if (flag == "--measurement_barrier") {
      if (value.back() == ',')
        throw std::invalid_argument("empty measurement barrier");
      std::stringstream in(value);
      std::string token;
      while (std::getline(in, token, ',')) {
        auto n = count(flag, token);
        if (n > std::numeric_limits<int>::max())
          throw std::invalid_argument("barrier exceeds range");
        option.measurement_barrier.push_back(static_cast<int>(n));
      }
    } else
      throw std::invalid_argument("unknown option: " + flag);
  }
  return option;
}

void cout_options(const Option &option) {
  std::cout.precision(16);
  std::cout << "# backend: " << option.backend << std::endl;
  std::cout << "# circuit: " << option.circuit << std::endl;
  std::cout << "# sparsepauli: " << option.sparsepauli << std::endl;
  std::cout << "# measurement barriers:";
  for (const auto step : option.measurement_barrier) {
    std::cout << " " << step;
  }
  std::cout << std::endl;
  std::cout << "# max_bp_iterations: " << option.max_bp_iterations << std::endl;
  std::cout << "# bp_tolerance: " << option.bp_tolerance << std::endl;
  std::cout << "# max_bond_dim: " << option.max_bond_dim << std::endl;
  std::cout << "# sv_min: " << option.sv_min << std::endl;
  std::cout << "# truncation_error: " << option.truncation_error << std::endl;
  if (!option.savename.empty()) {
    std::cout << "# savename: " << option.savename << std::endl;
  }
  if (!option.loadname.empty()) {
    std::cout << "# loadname: " << option.loadname << std::endl;
  }
  if (option.step_start != 0) {
    std::cout << "# step_start: " << option.step_start << std::endl;
  }
  if (option.step_end != 0) {
    std::cout << "# step_end: " << option.step_end << std::endl;
  }
}
