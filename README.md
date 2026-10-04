# Parallelized tensor network simulator with belief propagation

This repository provides a parallelized simulator for estimating expectation values of quantum circuits using belief propagation on tensor networks, built on top of the Tensor Computing API (TCAPI).

The library targets tensor-product-state representations—primarily MPS and PEPS—and provides utilities to construct circuit-induced factor graphs, run BP (loopy BP where appropriate), and aggregate observables efficiently across distributed computing resources.

**Key points**
- Focus: expectation-value evaluation of quantum circuits via tensor network simulation with belief propagation
- Built as an application-layer library on top of the Tensor Computing API (TCAPI)
- Tensor network formats: tensor product states such as MPS and PEPS
- Parallelism: thread-level and MPI-style process-level parallel execution
- Use cases: benchmark typical circuits/problems for paper-ready test calculations

> Note: This repository is intended as a foundation for test codes and benchmarking; the API may evolve as features are added.

## Tensor Computing Interface (TCI)

This project was originally developed in connection with the Tensor Computing Interface (TCI) project.

For details, see:

Sun, R.-Y., Shirakawa, T., Kohshiro, H., Sheng, D. N., Yunoki, S.  
*Tensor Computing Interface: An Application-Oriented, Lightweight Interface for Portable High-Performance Tensor Network Applications*  
https://arxiv.org/abs/2512.23917

The original TCI layer is being migrated to the public TCAPI specification.

### Current TCAPI integration

TNBP headers and sample programs use `tcapi/tcapi.h` and namespace `tcapi`.
The bundled CPU adapter is in `external/min-tci/include/tcapi`; its directory
name is retained. See [adapter documentation](external/min-tci/README.md) for
supported operations, binary I/O, diagnostics and limitations. The legacy
`external/min-tci/include/tci` is retained only for adapter compatibility checks;
TNBP does not maintain a second legacy implementation.

SVD/eigh outputs are now diagonal matrices. Tensor copies use `tcapi::copy`,
and linear combinations take `std::cref` inputs. `TensorProductState` construction
from tensors now takes `ctx` as the first argument; deep copying is explicit
with `state.copy(ctx)`, while ordinary copying is disabled and moves are allowed.

The validated backend is gqten on CPU. The historical Cytnx build branches and
external backend configurations are not validated for this TCAPI migration;
CUDA integration is a separate step. Existing CPU TPS stream format is retained.

## Requirements

To use the bundled TCAPI adapter, you need to install **GraceQ/tensor-ng-dev**.
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

### 1. Simulation of estimator for kicked ising model on a lattice (apps/kicked_ising_lattice)

The sample program under
```
/apps/kicked_ising_lattice
```
provides a simulator to estimate the expectation value after the kicked Ising Floquet dynamics on the lattice.
See the `/apps/kicked_ising_lattice/README.md` for more details.

### 2. Simulation of estimator for sparse Pauli operator after gate operations defined by qasm file (apps/estimator)

The apps program under
```
/apps/estimator/
```
provides a simulator for estimating the expectation value of sparse Pauli operator after gate operations defined by qasm file.
See the `/sample/estimator/README.md` for more details.

## Citation

If you use this code in your research, please cite:
```bibtex
@misc{Sun2025TCI,
  title        = {Tensor Computing Interface: An Application-Oriented, Lightweight Interface for Portable High-Performance Tensor Network Applications},
  author       = {Rong-Yang Sun and Tomonori Shirakawa and Hidehiko Kohshiro and D. N. Sheng and Seiji Yunoki},
  year         = {2025},
  eprint       = {2512.23917},
  archivePrefix= {arXiv},
  primaryClass = {quant-ph}
}
```

## Authors

- Tomonori Shirakawa
- Rongyang Sun
- Hidehiko Kohshiro

