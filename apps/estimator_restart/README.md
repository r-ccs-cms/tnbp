# Sample code for simulating the estimation of the expectation value of a sparse Pauli operator in the kicked-Ising model

## Overview

This repository provides a simplified simulator that performs two-dimensional tensor network simulations based on belief propagation along a given device topology. The simulator takes as input
1. a **QASM file** describing the quantum circuit, and
2. a **sparse Pauli operator file** for expectation value evaluation.
By combining these inputs, the simulator computes expectation values of quantum observables for models such as the transeverse and longitudinal field Ising Hamiltonian.

## Contents

- `estimator`
  The simulator reads the above QASM files and sparse Pauli operator data, and then executes belief-propagation-based two-dimensional tensor network simulations. Expectation values are computed as the final output.
  

## Installation

### 2. Estimator (c++ implementation)

The simulator is implemented on top of the Tensor Computing Interface (TCI),
which provides a unified API for tensor operations across different backends.

In order to run, both of the following components are required:

**(a) TCI interface** (`external/minitci`)
This repository includes `/external/minitci`, a header-only implementation of TCI that provides the interface layer required by the simulator.
It is included directly in this repository, so no additional setup is required.

**(b) Tensor backend** (`GraceQ/tensor-ng-dev`)
The actual tensor computations are carried out by GraceQ/tensor-ng-dev.
This library is included as a git submodule. To fetch it,
```
git submodule init
git submodule update
```
This will place the source under `external/tensor-ng-dev/`

**Installing GraceQ/tensor-ng-dev with the to-tci-release branch**
Go to to-tci-relase branch home directory, and type the following command
```
cmake -S . -B build \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
```
to setup the build, and type
```
cmake --build build -j
```
to build the hptt.

### 3. Building the Simulator

The simulator is **parallelized in real space** and therefore requires **MPI**.
Make sure you have an MPI environment available (e.g., OpenMPI, MPICH) and use an MPI C++ compiler such as mpicxx when building.

The build process is managed via the provided `Makefile`, which relies on the `Configuration` file for compiler and library paths.

By default, the Configuration file is set up for macOS environments where **LLVM (installed via Homebrew)** is used, including its OpenMP support. For example, the following flags are configured:
- `-fopenmp` and `-I/opt/homebrew/opt/llvm/include` (OpenMP include path)
- `-stdlib=libc++` (C++ standard library on macOS with LLVM)
- `-O3` (optimization)
The paths to the required libraries are already specified in `Configuration`:
- **TCI interface**: `external/minitci`
- **Tensor backend**: `external/tensor-ng-dev`
- **HPTT library**: bundled within `tensor-ng-dev/external/hptt`
- **TNBP headers**: relative path to this repository.
To build the simulator, simply run:
```
make
```
This will:
1. Use `mpicxx` as the c++ compiler
2. Include the header-only TCI implementation from `external/minitci`.
3. Link against the installed tensor backend (`tensor-ng-dev`), `LAPACK/BLAS`, and other system libraries specified in the configuration.
If you are on macOS with Homebrew-installed LLVM, no further changes are required.
For other environments (Linux clusters, alternative compilers, etc.), you may need to update the Configuration file to point to the correct OpenMP, LAPACK/BLAS, and MPI installations.

**Other environments**
For other systems (e.g., Linux clusters, Fujitsu A64FX, x86 HPC systems), you can adapt the build by editing the **Configuration** file.
We recommend keeping multiple configurations as commented-out blocks inside the file and switching them according to your environment.


## Usage

### Estimator
This program performs a real-space parallelized two-dimensional tensor network simulation based on belief propagation.
It takes as input a QASM circuit file and a sparse Pauli operator file, and computes expectation values of observables on a given device topology.
The simulation runs in parallel using MPI, so you should launch it with mpirun or mpiexec.

**Example**:
```
mpirun -np 4 ./estimator \
  --backend ibm_kobe \
  --circuit circuit.qasm \
  --sparse_pauli hamiltonian.dat \
  --max_bp_iterations 50 \
  --bp_tolerance 1.0e-8 \
  --max_bond_dim 100 \
  --sv_min 1.0e-8 \
  --truncation_error 1.0e-8
  --savename tnsdata
```

**Options**:
- `--backend <str>`:
  Target backend name (default: ibm_kobe).
- `--circuit <str>`:
  Input QASM file describing the circuit (default: circuit.qasm).
- `--sparse_pauli <str>`:
  Input file containing the sparse Pauli operator (default: sparsepauliop.txt).
- `--num_gates <list>`:
  Comma-separated list of the number of gates per layer (e.g., 4,4,4,4). If this list is not specified, the circuit automatically divides the Tensor Product Operator (TPO) into layers according to the barrier lines written in the QASM file.
- `--max_bp_iterations <int>`:
  Maximum number of belief propagation iterations (default: 50).
- `--bp_tolerance <float>`:
  Convergence tolerance for belief propagation (default: 1.0e-4).
- `--max_bond_dim <float>`:
  Maximum bond dimension 
- `--sv_min <float>`:
  Minimum cutoff of singular value to define the safe inverse.
- `--truncation_error <float>`:
  Target truncation error.
- `--savename <str>`:
  Filename to dump the tensor network state data
- `--loadname <str>`:
  Filename to load the tensor network state data for restart
