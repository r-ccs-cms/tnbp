# Parallelized tensor network simulator with belief propagation

This repository provides a parallelized simulator for estimating expectation values of quantum circuits
using **belief propagation** on tensor networks. The library targets tensor-product-state representations---primarily **MPS** and **PEPS**-- and exposes utilities to build circuit-induced factor graphs, run BP (loopy BP where appropriate), and aggregate observables efficiently across distributed resources.

**Key points**
- Focus: expectation-value evaluation of quantum circuit via tensor network simulation with belief propagation
- Tensor network formats: tensor product states such as MPS and PEPS
- Parallelism: thread-level and MPI-style process-level parallel executation
- Use cases: benchmark typical circuits/problems for paper-ready test calculations

> Note: This repository is intended as a foundation for test codes and benchmarking; the API may evolve as features are added.

## Tensor Computing Interface (TCI)

This repository provides a lightweight version of the **Tensor Computing Interface (TCI)**, derived from [GraceQ/tensor-dev-ng](https://github.com/GraceQ/tensor-dev-ng).  
The simplified TCI implementation is placed under:
```
/external/minitci
```

## Requirements

To use this TCI, you need to install **GraceQ/tensor-ng-dev**.  
For details on installation, please refer to the `README.md` of [tensor-ng-dev](https://github.com/gracequantum/tensor-ng-dev).

The repository is registered as a Git submodule and can be initialized as follows:

```bash
git submodule init
git submodule update
```
This will place the original `tensor-ng-dev` under
```
/external/tensor-ng-dev
```

## Sample Programs

### 1. Simulation of estimator for kicked ising model on a lattice (sample/kicked_ising_lattice)

The sample program under
```
/sample/kicked_ising_lattice
```
provides a simulator to estimate the expectation value after the kicked Ising Floquet dynamics on the lattice.
See the `/sample/kicked_ising_lattice/README.md` for more details.

### 2. Simulation of estimator for sparce pauli operator after gate operations defined by qasm file (sample/estimator)

The sample program under
```
/sample/estimator/
```
provides a simulator for estimating the expectation value of sparse Pauli operator after gate operations defined by qasm file.
See the `/sample/estimator/README.md` for more details.


## Authors

- Tomonori Shirakawa
- Rongyang Sun
- Hidehiko Kohshiro

