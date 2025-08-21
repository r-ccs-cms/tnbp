# Parallelized tensor network simulator with belief propagation

This repository provides a parallelized simulator for estimating expectation values of quantum circuits
using **belief propagation** on tensor networks. The library targets tensor-product-state representations---primarily **MPS** and **PEPS**-- and exposes utilities to build circuit-induced factor graphs, run BP (loopy BP where appropriate), and aggregate observables efficiently across distributed resources.

**Key points**
- Focus: expectation-value evaluation of quantum circuit via tensor network simulation with belief propagation
- Tensor network formats: tensor product states such as MPS and PEPS
- Parallelism: thread-level and MPI-style process-level parallel executation
- Use cases: benchmark typical circuits/problems for paper-ready test calculations

> Note: This repository is intended as a foundation for test codes and benchmarking; the API may evolve as features are added.


**Author**

- Tomonori Shirakawa
- Rongyang Sun
- Hidehiko Kohshiro

