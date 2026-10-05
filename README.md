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

The validated backends are gqten on CPU and tcapi-cuda for the estimator
applications and backend integration tests. See [CUDA setup](apps/README.cuda.md).
The historical Cytnx build branches are not validated for this TCAPI migration.
Existing CPU TPS stream format is retained.

## Requirements

### Backend repository access

As of October 2026, neither **tensor-ng-dev (gqten)** nor **tcapi-cuda** is
publicly available. Access permission to the selected backend repository is
required to clone its submodule and build with it. Listing these submodules in
TNBP does not grant access to their contents. Without permission, submodule
initialization will fail; a recursive clone attempts both restricted backends.

Initialize only the backend you intend to use, after obtaining access:

```sh
# CPU: gqten, used by the bundled adapter in external/min-tci
git submodule update --init --recursive external/tensor-ng-dev
# GPU: tcapi-cuda, pinned to the tested commit by this repository
git submodule update --init --recursive external/tcapi-cuda
```

For CPU installation, see [tensor-ng-dev](https://github.com/gracequantum/tensor-ng-dev).
For CUDA and cuTENSOR requirements and application build instructions, see
[CUDA setup](apps/README.cuda.md). A CUDA build does not require the gqten
submodule; a CPU build does not require the tcapi-cuda submodule.

## Sample Programs

- [Estimator](apps/estimator/README.md): apply a QASM circuit and evaluate local
  Pauli expectation values, using CPU or CUDA tensors.
- [Estimator with checkpoints](apps/estimator_restart/README.md): select
  measurement layers, repeat a circuit and save/load rank-local states.
- [Kicked Ising lattice example](apps/kicked_ising_lattice/README.md): simulate
  kicked Ising Floquet dynamics on a lattice.

For GPU setup and launch instructions for the two estimator applications, see
[the CUDA guide](apps/README.cuda.md). For small CPU/CUDA regression checks,
see [backend integration tests](test/backend_integration/README.md).

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

