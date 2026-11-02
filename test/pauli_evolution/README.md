# Pauli imaginary-time operator tests

`pauli_evolution` checks the tensor-independent route/layer plan and the TCAPI
exponential TPO builders. Configure through `test/` or this directory with the
same backend settings as `test/backend_integration`.

The independent reference applies Pauli strings to computational-basis bit
indices, enumerates the generated TPO virtual bonds on the host, and compares
complete matrices for up to five sites. Cases include Y signs, reversed and
unsorted edges, long-range strings, branching routes, cycles, isolated sites,
negative coefficients, zero/tiny/large time steps, and site-disjoint layers.
Weight-layout checks verify equal magnitudes at every non-identity Pauli site,
unit-weight I routers and branch points, a single negative sign at the root,
and weighting a Pauli branch point once rather than once per leg. These checks
run for complex float/double, alongside full-matrix references and local/full
tensor byte comparisons. A routed `1e-300` double-precision exponential checks
that sharing a tiny branch coefficient preserves its nonzero matrix element.
Commuting terms with shared support or intermediate identity routing sites
must be assigned to separate layers. Complex double and float
construction are exercised. Separate input and planning checks cover invalid
labels, coefficients, graphs and disconnected supports.

A dense Taylor exponential of a noncommuting two-site Hamiltonian checks that
halving the time step reduces first-order error by approximately two and
second-order error by approximately four. The identity energy shift is included
once per complete step. A two-site state test feeds a generated Y-containing
TPO to the existing `AbsorbTPO`; CPU MPI 1/2-rank tests cover local and distributed
site ownership. TPO generation itself performs no MPI communication.

```sh
cmake --build BUILD --target pauli_evolution
ctest --test-dir BUILD -R '^pauli_evolution_[12]_rank$' --output-on-failure
```

This suite validates operator construction and Trotter ordering, not BP ground
state convergence. CUDA builds use the existing test backend selection and are
launched explicitly inside a GPU allocation; each rank sees one GPU. Multiple
ranks may share that GPU in the shared-GPU configuration described below.

## CLI acceptance

`imaginary_time_cli` runs the application through its input and portable I/O
interfaces. It covers zero/plus/file-based states, staged schedules and both
Trotter orders, one-site/adjacent two-site energy including Y and identity terms,
measurement skipping and interval statistics, a thousand one-site steps,
state-only/warm-start/normal restart, changed MPI ownership, and invalid input or
corrupt/incomplete resume rejection. A native CPU record reader independently
contracts saved small states. Same-rank continuous/restarted states and
measurement-cadence variants must preserve V/E bytes.
Resume version 3 also records the operator-weight layout. Legacy version 2 and
different-layout metadata must not auto-continue; a valid legacy checkpoint used
as a new initial state must agree with state-only input. Missing layout metadata
in version 3 is rejected.

The repeated long-range example and branching/multiple-tree examples compare
saved states against independent exact bit-basis Trotter evolution. These are
small tree cases with BP between leaf-side batches, not a guarantee of general
exact BP compression. Metrics retain the maximum state distances. Tests cover
batch order, same-rank V/E restart identity, changed MPI ownership, connected
single-site normalization, and the versioned compression resume protocol.

## Batch compression API

`truncation_batch` tests two independent three-site operators on six sites:
`Z0 Y1 Z2` and `Y3 Z4 Y5`, plus routed two-site strings on the same components.
It groups their first edges into `{0,2}` and their second edges into `{1,3}`.
The caller's BP loop between batches gives small-system states matching an
independent bit-basis reference over ten steps. A frozen-environment control
records the existing approximation without asserting general exactness.

MPI 1/2/3/4-rank cases include remote endpoints, two independent rank pairs,
an empty rank and uninvolved ranks deliberately omitting a batch call. This
checks that the batch API introduces no global collective or barrier. Tests
also cover unselected V/E byte preservation, local selected-only results, empty
selection, invalid/duplicate/conflicting IDs, malformed local tensor counts, optional normalization accounting,
and equivalence of a no-BP batch sweep to the all-edge routine.

```sh
cmake --build BUILD --target truncation_batch
ctest --test-dir BUILD -R '^truncation_batch_[1234]_rank$' --output-on-failure
```

The numerical and endpoint-owner communication kernel is shared by both APIs.
The imaginary-time CLI calls the automatic leaf-side batches, refreshing BP
between them. Tests of all-edge equivalence with frozen environments still
exercise the unchanged all-edge API separately.

## Graph input compatibility

`graph_input` compiles the original two-argument graph reader as a function
pointer, verifies ascending endpoints without reordering edge IDs, checks
legacy sorting/deduplication, rejects malformed/overflowing endpoints and a
directory used as input, and accepts blank/comment-only rows.
CLI acceptance also covers that read failure and exact option-name validation.

## Leaf-side CLI connection review checks

The schedule contract is checked independently by removing each tree edge and
identifying its descendant component: every descendant must appear in an
earlier batch, each tree contributes at most one edge per batch, and selected
edges have no common endpoint. This includes all 24 edge-order permutations of
a five-site branching tree, independent trees, shared-site layer splitting and an
empty single-site schedule.

A zero BP-iteration budget exposes warnings before batches 0, 1 and 2 of the
branch example, checking that no final or measurement BP is inserted. Resume
coverage includes MPI2 continuous/restarted V/E byte identity, MPI2-to-MPI1
continuation against the exact state, and zero-evolution saves on changed owners.

Site-disjoint plan version 2 is also checked through the CLI: two strings with
disjoint Pauli support and disjoint edges sharing only the central identity
site are placed in separate layers. Exact evolution, same-rank V/E restart
identity and changed-rank continuation are checked. Previous plan versions are
rejected instead of silently changing the Trotter grouping on resume.

## Local TPO generation

Both term/layer APIs are checked against all-site generation byte-for-byte for
reordered subsets, a remote root, isolated sites, empty lists and zero/tiny/large
time steps. Invalid or duplicate targets are rejected. Whole-operator scales
are identical for every selection; complex float/double are covered.
Local `AttachTPO` and global `AbsorbTPO` must produce identical V/E records.
The repeated batch evolution tests use local generation on MPI1/2/3/4 ranks,
including an empty rank and independent owner pairs. CLI exact evolution and
restart acceptance exercise the local path throughout.

Local-layout review adds an unsorted branching graph with a remote root,
noncontiguous owners, reversed local site/V order and reversed local edge/E
order. Two successive attachments must match global absorption byte-for-byte;
active message dimensions grow from 1 to 2 to 4 while the unused boundary edge
stays at dimension 1. Both directed messages are checked against normalized
identity matrices, and an isolated local site is included. CLI checkpoints also
check whole-operator log scale against an independent coefficient sum before
and after continuation on different rank counts.

## CUDA validation

On November 1, 2026, small-system tests passed on ROQUO NVIDIA GB200 GPUs
with tcapi-cuda `f32c87d`: one shared GPU with MPI 1/2/3/4 ranks for batch
compression and MPI 1/2 ranks for operator tests, followed by same-node
GPU-per-rank execution with MPI 1/2 ranks. Both configurations passed CLI
acceptance (74 application runs and 2239 checks each), covering initial states,
time schedules, measurement/statistics and save/resume. The GPU acceptance
runner used a CUDA-native checkpoint reader in place of the CPU-native reader
in `cli.py`; reference calculations and numerical tolerances were unchanged.
The largest exact-state distance across CLI cases was below `2e-15`. A checkpoint
saved in the shared-GPU job and repacked in the two-GPU job was byte-identical.

These results supplement CPU validation; they do not establish general BP
convergence, multi-node correctness or performance. CPU/CUDA native tensor
checkpoint interchange is not supported. CPU-only CTest registrations remain
disabled for CUDA builds.

Support-balanced weights were revalidated on 2026-11-02 at commit `cdda0c2`,
using the same CUDA backend `f32c87d`. Both shared-GPU and two-GPU runs passed
complex float/double operator construction on MPI 1/2 and 78 CLI runs with
2,375 checks each. Batch tests passed shared-GPU MPI 1/2/3/4 and distinct-GPU
MPI 1/2. A shared-GPU checkpoint repacked on two distinct GPUs was byte-identical.
The committed CLI acceptance script was used with only its native-record decoder
adapted for GPU records. This supersedes the earlier 74-run count for the new
weight layout; the validation scope limitations above still apply.
