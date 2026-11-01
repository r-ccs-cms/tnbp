# Pauli imaginary-time operator tests

`pauli_evolution` checks the tensor-independent route/layer plan and the TCAPI
exponential TPO builders. Configure through `test/` or this directory with the
same backend settings as `test/backend_integration`.

The independent reference applies Pauli strings to computational-basis bit
indices, enumerates the generated TPO virtual bonds on the host, and compares
complete matrices for up to five sites. Cases include Y signs, reversed and
unsorted edges, long-range strings, branching routes, cycles, isolated sites,
negative coefficients, zero/tiny/large time steps, and site-disjoint layers.
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
state convergence. CUDA builds use the existing test backend selection and must
be launched explicitly with one visible GPU per rank. CPU/CUDA agreement is not
implied by CPU-only results.

## CLI acceptance

`imaginary_time_cli` runs the application through its input and portable I/O
interfaces. It covers zero/plus/file-based states, staged schedules and both
Trotter orders, one-site/adjacent two-site energy including Y and identity terms,
measurement skipping and interval statistics, a thousand one-site steps,
state-only/warm-start/normal restart, changed MPI ownership, and invalid input or
corrupt/incomplete resume rejection. A native CPU record reader independently
contracts saved small states. Same-rank continuous/restarted states and
measurement-cadence variants must preserve V/E bytes.

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
