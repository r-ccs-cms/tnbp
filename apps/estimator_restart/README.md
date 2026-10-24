# Quantum-circuit estimator with checkpoints

`estimator_restart` adds rank-local checkpoint loading/saving and measurements
at selected TPO layers. Each step applies the entire supplied QASM circuit;
multiple steps repeat that circuit. The executable in this directory is named
`estimator`, as in the ordinary estimator directory.

## Portable checkpoint input and output

See [portable TPS checkpoints](../../docs/portable_checkpoints.md) for independent
site/message files shared by estimator and sampler, including restart with a
new MPI process count. This is a separate format from the legacy rank-local
options described below. No external library dependency is added.

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
once. These are repetition indices, not QASM instruction offsets. Legacy rank-local loading requires an explicit next step. Portable loading
uses the app-owned optional `estimator.resume` file to infer the saved next step
for identical QASM contents, or starts at zero for a different circuit while
retaining the saved state and graph. Without that file, state-only input starts
at zero. The library manifest contains no circuit or step metadata. An explicit
`--step_start` must match the inferred start. See the portable checkpoint guide
for state preparation followed by a different evolution circuit.

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
| `--portable_save` | Unset | Save an immutable portable TPS directory after the complete cycles. |
| `--portable_load` | Unset | Load a portable TPS and validate estimator resume metadata on the new MPI partition. |
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

## BP convergence

BP stops when the mean directed-edge Frobenius distance between the previous
and updated messages is below `--bp_tolerance`. Messages have unit Frobenius norm and positive-real trace at initialization,
attach, truncation, restart loading and BP update boundaries. The check only
subtracts successive messages and sums their Frobenius distances. The
updated messages are committed even on the terminating iteration. Each
directed edge contributes on its source site's owner and the scalar residual
is combined with MPI Allreduce.

This differs from the legacy fixed-point check, which constructs another
prediction from the updated state. Iteration counts and final values can
therefore differ at the same tolerance. The library keeps
`BeliefPropagationCondition` for that legacy check; these app loops use
`BeliefPropagationDifference`. No site-sized prefix cache is retained.

## Site ownership

Use `--ownerfile FILE` for `site_id mpi_rank` rows; omitted selects the
automatic contiguous partition. The selected assignment is printed once as
`# site ownership: ...` and `# MPI rank N: ...` without timestamps. See
[site ownership](../../docs/site_ownership.md) for input validation, direct
initialization/loading and legacy checkpoint restrictions.

## Portable publication and gate topology

The estimator app stages the state and its resume file together under
`<destination>.estimator.pending/state` and publishes the complete directory
with one rename. A failure before publication leaves no new final checkpoint; its staging
files remain for inspection. A legacy leftover `estimator.resume.pending`
is an error for estimator input. Library-only state input remains supported.
TPO construction rejects a gate whose operands cannot be connected on the saved
graph; disconnected graphs themselves and gates within their components remain
valid. See [portable checkpoints](../../docs/portable_checkpoints.md).

This publication path and unrepresentable-gate handling passed 19 dedicated
regression cases on CPU and on each of ROQUO's shared-one-GPU and two-GPU
configurations. Tests observe a save between state completion and app metadata
creation, and inject failures before metadata, during metadata writing and
before final publication. Valid component-local gates and connected-graph
routing remain supported.
