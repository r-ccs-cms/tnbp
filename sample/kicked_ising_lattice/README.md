# Sample code for simulating the estimate of the expectation value of sparse Pauli operator in kicked Ising dynamics

## Overview

This sample proides a simplified simulator that peforms one and two-dimensional tensor network simulations based on belief propagation. Namely,
$\langle 0 \vert \hat{U}^{\dagger} \hat{O} \hat{U} \vert 0 \rangle$ where $\hat{U} = ( \prod_{i} \exp ( -i \pi h_z \tau \hat{Z}_i ) \prod_{\langle i,j \rangle} \exp ( - i \pi J_z \tau \hat{Z}_i \hat{Z}_j ) \prod_{i} \exp ( - i \pi h_x \tau \hat{X}_i ) )^s \vert 0 \rangle$, $\hat{O} = \{ \{ \hat{X}_i \}_i, \{ \hat{Z}_i \}_i, \{ \hat{Z}_i \hat{Z}_j \}_{\langle i,j \rangle}$.

## Installation

### 1. Simulator (c++ implementation)

This simulator is implemented on top of the Tensor Computing Interface (TCI), which provides a unified API for tensor operations across different backends.

In order to run, both of the following components are required:

**(a) TCI interface** (`external/minitci`)
This repository includes `/external/minitci`, a header-only implementation of TCI that provides the interface layer required by the simulator.
It is included directly in this repository, so no additional setup is required.

**(b) Tensor backend** (`GraceQ/tensor-ng-dev`)
The actual tensor computations are carried out by `GraceQ/tensor-ng-dev`.
This library is included as a git submodule. To fetch it,
```
git submodule init
git submodule update
```
This will place the source under `external/tensor-ng-dev/`

**Installing GraceQ/tensor-ng-dev**
Follow the official installation guide.
A simplified installation command is:
```
make install \
     CXX="your_preferred_compiler" \
     CMDLINE_CXX_FLAGS="include paths for OpenMP or other dependencies" \
     CMDLINE_LINK_FLAGS="linker flags for LAPACK, BLAS, etc."
```
For example, in a macOS environment, when using `clang++` installed via `brew install llvm`, we add the following options:
```
make install \
     CXX="/opt/homebrew/opt/llvm/bin/clang++" \
     CMDLINE_CXX_FLAGS="-std=c++20 -stdlib=libc++ -fopenmp -I/opt/homebrew/opt/llvm/include -O3" \
     CMDLINE_LINK_FLAGS="-L/opt/homebrew/opt/openblas/lib -llapack -lblas"
```

### 2. Building the simulator

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

This program performs a real-space parallelized two-dimensional tensor network simulation based on belief propagation.
It takes as input a QASM circuit file and a sparse Pauli operator file, and computes expectation values of observables on a given device topology.
The simulation runs in parallel using MPI, so you should launch it with mpirun or mpiexec.

**Example**:
```
mpirun -np 4 ./kicked_ising_exvalue \
  --lattice oned \
  --L 8 \
  --Jz 0.5 \
  --hz 0.0 \
  --hx 0.9 \
  --dt 0.5 \
  --num_steps 50 \
  --max_bp_iteration 50 \
  --bp_tolerance 1.0e-8 \
  --max_bond_dim 100 \
  --sv_min 1.0e-6 \
  --truncation_error 1.0e-6
```

**Options**:
- `--latice <str>`: Target lattice. Currently, `oned` and `honeycomb` are available.
- `--L <int>`: The size of system. For `oned`, it is the system size. For `honeycomb`, the system size becomes $2 \times L^2$.
- `--Jz <float>`: The coupling strength for $\hat{Z}_i \hat{Z}_j$ term.
- `--hx <float>`: The strength of external field along x-direction.
- `--hz <float>`: The strength of external field along z-direction.
- `--dt <float>`: Time slice.
- `--num_steps <int>`: Number of steps of Floquet cycles.
- `--max_bp_iteration <int>`: Maximum number of iteration for belief propagation procedure.
- `--bp_tolerance <float>`: Target tolerance for iteration for belief propagation procedure.
- `--max_bond_dim <int>`: Maximum bond dimension for virtual bond.
- `--sv_min <float>`: Minimum cutoff in singular value.
- `--truncation_error <float>`: Target truncation error.
