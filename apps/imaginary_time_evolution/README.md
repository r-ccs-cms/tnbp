# Imaginary-time evolution with belief propagation

`imaginary_time_evolution` applies a Pauli-string Hamiltonian's imaginary-time
Trotter operators to a tensor-network state, using the estimator BP/truncation
pipeline. Second-order symmetric Trotter splitting is the default. A layer
contains terms whose complete routing trees share no sites, including
intermediate identity sites. Arbitrary
many-site and long-range Pauli strings can be evolved when their support is
connected in the supplied graph. Measurement supports one-site and graph-edge
two-site terms; general-string measurement is not implemented.

## Build

Use the CPU/CUDA dependencies and configuration described in
[estimator](../estimator/README.md) and the [CUDA guide](../README.cuda.md).
Build with `make CONFIGURATION=/path/to/Configuration`, or configure this
application directory with the same CMake backend parameters used by the tests.
The common `test/` build also provides the `imaginary_time_evolution` target.

## New calculation

The Hamiltonian text file uses `real [imag] PauliString` rows, for example:

```text
0.23 0 IY
-0.37 0 XI
0.70 0 ZZ
0.11 0 II
```

All strings have the same nonzero length N; the rightmost character is site 0.
Equal labels are merged and their resulting coefficients must be exactly real.
Nonfinite coefficients and unsupported labels are rejected. Identity terms are
included in energy and tracked as a scalar factor, not applied to the state.

The graph file has one `site_a site_b` row per undirected edge. Endpoints are normalized to ascending site order; the edge row order is
preserved. See the [shared graph conventions](../../docs/graph_tensor_conventions.md). Duplicate edges, self loops and endpoints
outside `0..N-1` are rejected. All N sites are retained, including isolated sites.
An empty graph file supports independent one-site problems. Multi-site terms
spanning disconnected components are rejected.

```sh
mpirun -np 1 ./imaginary_time_evolution \
  --inputedge edges.txt --sparse_pauli H.txt \
  --initial_state plus --time_step 0.05 --num_steps 100 \
  --measure_every 10 --portable_save checkpoint-100
```

The initial state defaults to all zero. `--initial_state plus` selects the
uniform product state. Alternatively, `--initial_statefile FILE` supplies exactly
one `site_id bit` row for every site, with bit 0 or 1; row order is arbitrary and
`#` comments are allowed. These initial-state options are mutually exclusive.
Neither fixed product state guarantees overlap with the ground space.

For several time-step stages, replace `--time_step`/`--num_steps` with
`--schedule schedule.txt`:

```text
# time_step num_steps
0.10 20
0.03 40
0.01 100
```

This schedule has 160 complete Trotter steps and total imaginary time 4.2.
Each step includes the full forward/reverse layer sequence. After each layer
action, BP converges and the used edges are compressed in leaf-side batches
with BP between batches. There is at most one edge from each term tree per
batch, with no shared endpoints across trees. Ready edges use global edge ID
order; conflicting tree heads wait. There is no BP after the final batch.
See [automatic batch scheduling](../../docs/truncation_batches.md#automatic-leaf-side-schedules-for-pauli-layers). Half-step boundaries are not merged across
steps. There is no energy-based automatic stopping rule.

The graph, plan and operator scalar are calculated on every rank. TPO tensors
are generated and cached only for that rank's local sites, in the same order
as its state tensors, and applied through `AttachTPO`. `log_scale` is a common
whole-operator value, not a rank-local contribution to sum. Local generation
does not change native tensor records or the plan/compression protocol identifiers.

The TPO's branch coefficients are shared equally over the non-identity Pauli
sites. Intermediate I sites, including routing branch points, carry no weight;
the sign of the Pauli-branch coefficient appears once at the term root. This
`pauli-support-equal-v1` layout leaves the exact exponential and `log_scale`
unchanged and preserves the existing cutoff and parallel compression routines.

## Options and output

| Option | Default | Meaning |
|---|---|---|
| `--inputedge FILE` | Required for fresh state | Ordered graph edges; optional with portable input and must match if given. |
| `--sparse_pauli FILE` | Required for new calculation | Hamiltonian text input. |
| `--initial_state zero\|plus` | `zero` | Fresh product-state choice. |
| `--initial_statefile FILE` | None | Per-site computational basis bits. |
| `--time_step DT --num_steps N` | None | One schedule stage, both strictly positive. |
| `--schedule FILE` | None | Stage rows; exclusive with the preceding pair. |
| `--trotter_order 1\|2` | `2` | First or second order. |
| `--max_bp_iterations N` | `50` | Maximum updates per layer action; zero is allowed and reports nonconvergence. |
| `--bp_tolerance EPS` | `1e-8` | Mean directed-message change threshold. |
| `--max_bond_dim D` | `100` | Maximum retained bond dimension. |
| `--sv_min EPS` | `1e-8` | Existing truncation/root cutoff. |
| `--truncation_error EPS` | `1e-8` | Existing SVD truncation target. |
| `--measure_every K` | `1` | Measure at absolute completed steps divisible by K; K must be positive. |
| `--step_end K` | Schedule endpoint | Exclusive absolute endpoint, not an additional number of steps. |
| `--portable_save DIR` | None | Save a new checkpoint directory at normal completion. |
| `--portable_load DIR` | None | Load V/E and the saved ordered graph. |
| `--restart_mode auto\|initial` | `auto` | Continue this app's saved calculation or explicitly use its state as a new initial state. |
| `--ownerfile FILE` | Balanced contiguous | Existing `site_id rank` ownership layout. |
| `-h`, `--help` | — | Show help and exit. |

Normal completion always measures the endpoint once, even if it is not a
multiple of K. Measurement does not run additional BP and does not modify V/E.
A summary reports energy, imaginary time, current bond-dimension min/max, and
truncation-error min/max over all compression operations since the preceding
summary. The window restarts at the loaded step on restart. No edges yield
`N/A` bond dimensions; no compression operations yield `N/A` truncation errors. The reported truncation errors
are the existing per-edge SVD diagnostics, not global state-error bounds. In
particular, inverse-square-root cutoffs introduce support projections that are
not counted in those SVD tails. Combining projections from the same pre-update
environments can change the state even on a tree, despite converged BP and zero
SVD discarded weight. The CLI now refreshes BP between leaf-side batches, while preserving the
existing per-edge numerical kernel and cutoff. This improves the tested tree
examples but does not guarantee exactness for arbitrary states or loop graphs.
The original all-edge API remains unchanged. Simply lowering the cutoff need
not resolve the approximation.

For the distinction between evolution and measurable energy terms, see
[measurement limits](#measurement-limits).
Input information, including the compression protocol and each layer's batch
edge IDs, starts with `# <timestamp>`; results, warnings and controlled
library diagnostics also carry a timestamp. Normal operation does not emit
per-edge or per-BP-iteration progress. Explicit backend verbose modes and MPI
runtime diagnostics retain their backend/launcher behavior.

## Measurement limits

Leaf-side batch compression and BP do not extend the measurement API. The
complete Hamiltonian, including long-range and many-site terms, is still used
for imaginary-time evolution; unsupported measurement terms are skipped only
when forming the reported energy.

| Normalized Hamiltonian term | Evolution | Energy measurement |
|---|---|---|
| Identity | Tracked as a scalar | Added as a constant shift |
| One non-identity site | Supported | Supported |
| Two non-identity sites joined by a graph edge | Supported | Supported |
| Two non-identity sites not joined by a graph edge | Supported if connected by a route | Skipped |
| Three or more non-identity sites, including branching support | Supported if connected by a tree | Skipped |

Adjacency means an edge in the supplied graph, not consecutive site numbers.
For example, on `0--1--2`, a term acting on sites 0 and 2 is evolved through
both edges but is not measured. A three-site term is also not measured, even
when its entire route is a simple path or tree.

When any term cannot be measured, summaries print `energy=N/A`, the number of
skipped terms as `unsupported_terms`, and `partial_energy` for the sum of
supported terms plus the identity shift. This partial sum is not the full
Hamiltonian energy and must not be used as its ground-state convergence
criterion. If all non-identity terms are unsupported, the partial sum contains
only the identity shift (zero when absent).

The evolution, bond-dimension min/max and truncation-error min/max continue
normally. `--measure_every` controls measurement/summary frequency, with a final
summary at normal completion. Measurement uses the current messages without
additional BP. General Pauli-string energy measurement is not implemented.

## Save and resume

```sh
# Stop at a complete step, inside an otherwise longer schedule.
mpirun -np 1 ./imaginary_time_evolution \
  --inputedge edges.txt --sparse_pauli H.txt --schedule schedule.txt \
  --step_end 60 --portable_save checkpoint-60

# Restore schedule, Hamiltonian and BP/truncation conditions, on new owners.
mpirun -np 2 ./imaginary_time_evolution \
  --portable_load checkpoint-60 --step_end 160 \
  --portable_save checkpoint-160
```

Resume metadata version 3 records the `pauli-support-equal-v1` operator weights
and the unchanged `leaf-batches-v1` compression protocol. Plan version 2 still
requires site-disjoint trees. Version 2 resume metadata describes the former
root-only weights: automatic continuation is rejected rather than silently
changing the numerical trajectory. A valid version 2 checkpoint can still be
used with `--restart_mode initial` and an explicitly supplied Hamiltonian and
schedule. Earlier draft metadata or plan versions remain rejected. Generic
portable states without app metadata remain valid initial states; the native
V/E format is unchanged.

A valid `imaginary_time.resume` restores the plan identity, complete schedule,
calculation settings and next step. Explicit settings must match. MPI count,
owner assignment, measurement frequency, endpoint and output location may
change. Portable input requires at least one site per reader rank and a matching
backend/scalar codec. CPU/CUDA checkpoint interchange is not supported.

To deliberately change Hamiltonian, schedule or numerical settings, pass
`--restart_mode initial` and specify a new Hamiltonian and schedule. V/E are
retained, but step, imaginary time and the tracked operator scalar reset to zero.
A generic portable state without this app's resume file is also a new initial
state; metadata belonging to other apps is ignored. Corrupt/incomplete metadata
belonging to this app is rejected even in `initial` mode.

State and resume data are staged together, then published by directory rename.
The destination must not already exist and its parent must exist. No periodic
save, signal-time save or mid-step resume is provided. `step_end` equal to the
loaded next step measures/saves without evolution. The tracked operator log
scale accounts for deliberately removed exponential coefficients; it does not
include every state normalization and is not a partition function or state norm.

See [operator API](../../docs/pauli_evolution.md),
[portable format](../../docs/portable_checkpoints.md), and
[acceptance tests](../../test/pauli_evolution/README.md).

## Validation scope

Small-system CPU and ROQUO NVIDIA GB200 tests cover initial states, schedules,
measurement/statistics, leaf-side batch compression with intervening BP, and
save/resume. CUDA validation includes shared-GPU MPI 1/2/3/4-rank batch tests
and same-node GPU-per-rank MPI 1/2-rank tests; CLI acceptance passed 74 runs in
each GPU configuration. See [test coverage](../../test/pauli_evolution/README.md#cuda-validation)
for the backend revision and checkpoint-reader setup. These tests do not
establish general BP convergence or multi-node/performance guarantees.

The support-balanced weight layout was revalidated on GB200 with CUDA backend
`f32c87d`: operator construction for complex float/double on MPI 1/2, batch
compression with shared-GPU MPI 1/2/3/4 and two-GPU MPI 1/2, and 78 CLI
acceptance runs in each GPU configuration. Reading and repacking a shared-GPU
checkpoint on two distinct GPUs also preserved all checkpoint bytes.
