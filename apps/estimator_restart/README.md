# Quantum-circuit estimator with checkpoints

`estimator_restart` adds rank-local checkpoint loading/saving and measurements
at selected TPO layers. Each step applies the entire supplied QASM circuit;
multiple steps repeat that circuit. The executable in this directory is named
`estimator`, as in the ordinary estimator directory.

## Build

For CPU dependencies and configuration, follow the
[estimator build instructions](../estimator/README.md#build), then build from
`apps/estimator_restart`:

```sh
make
# Or use a configuration adapted to your installation:
make CONFIGURATION=/path/to/Configuration.cpu
```

For GPU builds and one-GPU-per-rank launching, follow the shared
[CUDA guide](../README.cuda.md). Both gqten and tcapi-cuda currently require
repository access permission; see [backend access](../../README.md#backend-repository-access).

## Circuit steps, layers and measurements

QASM barriers divide the circuit into TPO layers. `--measurement_barrier`
selects zero-based layer indices; measurements occur after applying BP and
truncation for those layers, on every executed step. If it is omitted, no
expectation values are printed.

`--step_start` is inclusive and `--step_end` is exclusive. With both options
omitted (or zero), the effective range is `[0, 1)`, so the circuit is applied
once. These are repetition indices, not QASM instruction offsets. Loading a
checkpoint does not infer the next step number: specify it explicitly.

Unlike the ordinary estimator, this application does not implement
`--num_gates`, `--do_opt_tpo` or `--eps_opt_tpo`.

## Example: save and continue

From `apps/estimator_restart`, the following uses the bundled kicked Ising
inputs from the neighboring estimator directory. That input has four TPO
layers, so index `3` selects the final layer.

```sh
# Apply the circuit once and save the final state.
mpirun -np 2 ./estimator \
  --backend ibm_kobe \
  --circuit ../estimator/kicked_ising.qasm \
  --sparse_pauli ../estimator/H_ising_kobe.txt \
  --measurement_barrier 3 \
  --step_end 1 \
  --savename checkpoint-

# Load that state, apply the same circuit once more and save to a new prefix.
mpirun -np 2 ./estimator \
  --backend ibm_kobe \
  --circuit ../estimator/kicked_ising.qasm \
  --sparse_pauli ../estimator/H_ising_kobe.txt \
  --measurement_barrier 3 \
  --step_start 1 --step_end 2 \
  --loadname checkpoint- \
  --savename continued-
```

This illustrates the continuation interface; the kicked Ising continuation
above has not been part of the current validation. For GPU execution, use the
launcher and binding described in the [CUDA guide](../README.cuda.md#gpu-assignment).

Checkpoints are saved once, after all requested steps. Each rank writes its
own file: for example, prefix `checkpoint-` produces `checkpoint-000000.dat`,
`checkpoint-000001.dat`, and so on. Provide an existing output directory and
use a new prefix to preserve previous files. Restart with the same tensor
backend, precision, MPI rank layout and graph/circuit setup. Cross-backend
checkpoint conversion is not supported. Keep the original input files and
step information alongside the checkpoints.

## Options

Every option takes a value. Supply valid nonnegative step/iteration counts and
valid layer indices; the CLI does not provide comprehensive validation.

| Option | Default | Meaning |
|---|---|---|
| `--backend` | `default` | Graph topology: QASM-derived edges (`default`) or built-in `ibm_kobe`. Not the tensor backend. |
| `--circuit` | `circuit.qasm` | Input QASM path; required even when loading a checkpoint. |
| `--sparse_pauli` | `sparsepauliop.txt` | Input Pauli path; read even if no measurement layers are selected. |
| `--measurement_barrier` | Unset | Comma-separated zero-based TPO layer indices, e.g. `0,3`. Unset means no measurements. |
| `--savename` | Unset | Prefix for checkpoint files saved at the end. |
| `--loadname` | Unset | Prefix for checkpoint files to load; otherwise initialize a product state. |
| `--step_start` | `0` | First circuit repetition index, inclusive. |
| `--step_end` | `0` (effective `1`) | End repetition index, exclusive; `0` selects the default `1`. |
| `--max_bp_iterations` | `50` | Maximum BP iterations per TPO layer. |
| `--bp_tolerance` | `1e-8` | BP convergence threshold. |
| `--max_bond_dim` | `100` | Maximum retained bond dimension (integer). |
| `--sv_min` | `1e-8` | Singular-value cutoff used during truncation. |
| `--truncation_error` | `1e-8` | Target truncation error. |

Input formats and observable limitations are the same as the ordinary
[estimator](../estimator/README.md#inputs-and-output). Output adds the step
and layer index to each expectation value. Values are individual Pauli
expectations: input coefficients are not applied, and terms are not summed.

## Validation

A two-qubit Bell case passed on CPU with two MPI ranks and on ROQUO with
one/two GPUs. Loading and saving without further circuit application preserved
checkpoint bytes for each rank. Continued nontrivial evolution, arbitrary
restart options and cross-backend checkpoint loading have not been validated
in this integration. See the [CUDA validation summary](../README.cuda.md#validation).
