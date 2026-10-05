# Quantum-circuit expectation estimator

`estimator` applies a QASM circuit to an initial tensor product state, using
belief propagation and bond truncation, then evaluates local Pauli observables.
It supports the gqten CPU backend and the tcapi-cuda GPU backend through TCAPI.

## Build

Both backend repositories currently require access permission. See
[backend repository access](../../README.md#backend-repository-access).
Python and Qiskit are only needed to generate new inputs, not to run the
C++ application with the bundled QASM and Pauli files.

### CPU / gqten

From the TNBP root, initialize the pinned CPU dependency:

```sh
git submodule update --init --recursive external/tensor-ng-dev
```

Build gqten's HPTT dependency following its own instructions. The application
uses the bundled header-only TCAPI adapter in `external/min-tci/include`;
there is no separate TNBP library to install.

A C++17 compiler, MPI, OpenMP, BLAS/LAPACK and HPTT are required. From
`apps/estimator`, adapt `Configuration` to your installation and run:

```sh
make
```

The supplied `Configuration` is a macOS/Homebrew example, not a portable
installation detector. Check the MPI wrapper's compiler, OpenMP flags,
include directories and library paths. In particular, HPTT is under
`Source/gqten/ext/hptt` in the backend source; the library location depends
on how it was built. A separate configuration file can be selected with:

```sh
make CONFIGURATION=/path/to/Configuration.cpu
```

### GPU / tcapi-cuda

Follow the shared [CUDA build and launch instructions](../README.cuda.md).
They cover the pinned submodule, CUDA/cuTENSOR paths and one visible GPU per
MPI rank. The tensor backend is selected at build time; the CLI `--backend`
option selects graph topology.

## Run the bundled example

From `apps/estimator`, after building the desired executable:

```sh
# CPU, one MPI rank; uses the existing input and calculation settings
sh run.sh
```

For multiple CPU ranks, use the same inputs and options with a different count:

```sh
mpirun -np 2 ./estimator \
  --backend ibm_kobe \
  --circuit kicked_ising.qasm \
  --sparse_pauli H_ising_kobe.txt \
  --max_bp_iterations 50 \
  --sv_min 1.0e-4 \
  --truncation_error 1.0e-6
```

For GPU execution, replace the MPI launcher with the site-specific GPU binding
shown in the [CUDA guide](../README.cuda.md#gpu-assignment).

## Inputs and output

- The QASM file supplies the circuit. By default, barriers divide it into TPO
  layers; `--num_gates` can instead specify gate counts per layer.
- The sparse Pauli file supports lines such as `1.0 0.0 ZIZI`, meaning real
  coefficient, imaginary coefficient and Pauli label. The rightmost character
  represents qubit 0.
- The current measurement path selects one-site terms and two-site terms on
  graph edges. Use these local observables; general many-site and identity-only
  terms are not supported by this application path.
- Output contains BP residuals, truncation errors, retained bond dimensions and
  individual Pauli expectation values. The current tensor conversion does not
  apply the input coefficients, and the application does not sum a Hamiltonian
  expectation. Apply coefficients and sum the terms separately when needed.

## Options

Every option below takes a value. Supply valid nonnegative iteration counts
and layer specifications; the CLI does not provide comprehensive validation.

| Option | Default | Meaning |
|---|---|---|
| `--backend` | `default` | Graph topology: `default` derives edges from QASM; `ibm_kobe` uses the built-in graph. |
| `--circuit` | `circuit.qasm` | Input QASM path. |
| `--sparse_pauli` | `sparsepauliop.txt` | Input Pauli file path. |
| `--num_gates` | Unset | Comma-separated gate counts per layer, e.g. `4,4`; otherwise use QASM barriers. |
| `--max_bp_iterations` | `50` | Maximum BP iterations per TPO layer. |
| `--bp_tolerance` | `1e-8` | BP convergence threshold. |
| `--max_bond_dim` | `100` | Maximum retained bond dimension (integer). |
| `--sv_min` | `1e-8` | Singular-value cutoff used during truncation. |
| `--truncation_error` | `1e-8` | Target truncation error. |
| `--do_opt_tpo` | `0` | Nonzero enables TPO optimization by SVD. |
| `--eps_opt_tpo` | `1e-8` | Tolerance for optional TPO optimization. |

## Generate new inputs (optional)

The Python utilities require Qiskit. Querying an IBM backend also requires
`qiskit_ibm_runtime` and configured service credentials. Choose a backend
available to your account; the bundled inputs require no service access.

```sh
pip install qiskit qiskit_ibm_runtime
python kicked_ising_qasm.py --backend ibm_kobe --steps 1 --output circuit.qasm
python dump_sparse_pauli.py --backend ibm_kobe --out observables.txt
```

`kicked_ising_qasm.py` accepts `--steps` (default `1`), `--hx` (`0.5`),
`--hz` (`0.0`), `--jz` (`0.5`), `--backend` (`ibm_kobe`) and `--output`
(`kicked_ising.qasm`). The rotation angles are these coefficients times pi.

`dump_sparse_pauli.py` accepts `--hx`, `--hz`, `--jz` (each default `1.0`),
`--out` (default `ising_topology.txt`), and optional `--ibm-instance`.
For an offline topology, supply `--n-qubits` and `--edge-file` (lines `i j`),
or `--n-qubits` and `--fully-connected`, instead of `--backend`. Keep the
circuit and observable qubit numbering/topology consistent.

## Validation

The bundled 156-qubit kicked Ising case completed on CPU with 1/2/4 MPI ranks
and on ROQUO with 1/2 GPUs. All 488 expectation values agreed within `8e-14`
with the CPU one-rank result. This is an execution and backend-consistency
check, not an independent exact solution or a performance benchmark.
See also the [backend integration tests](../../test/backend_integration/README.md).
