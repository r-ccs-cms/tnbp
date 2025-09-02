# Parallelized tensor network simulator with belief propagation

This repository provides a parallelized simulator for estimating expectation values of quantum circuits
using **belief propagation** on tensor networks. The library targets tensor-product-state representations---primarily **MPS** and **PEPS**-- and exposes utilities to build circuit-induced factor graphs, run BP (loopy BP where appropriate), and aggregate observables efficiently across distributed resources.

**Key points**
- Focus: expectation-value evaluation of quantum circuit via tensor network simulation with belief propagation
- Tensor network formats: tensor product states such as MPS and PEPS
- Parallelism: thread-level and MPI-style process-level parallel executation
- Use cases: benchmark typical circuits/problems for paper-ready test calculations

> Note: This repository is intended as a foundation for test codes and benchmarking; the API may evolve as features are added.

## Requirements

- **External tensor library (mandatory):** This project depends on **TCI (tensor computing interface)** as the core tensor computation backend.
  Please install and make it discoverable (e.g., via `CMAKE_PREFIX_PATH` or environment variables) before building this repository.

## Getting Started

```bash
# assuming the external tensor library is installed and discoverable
git clone https://github.com/r-ccs-cms/tnbp.git
cd tnbp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=/path/to/gqten/tci
cmake --build build -j
```

## Authors

- Tomonori Shirakawa
- Rongyang Sun
- Hidehiko Kohshiro

