# MCMC sampling tests

Build with the CPU configuration of the [common test entry point](../README.md#build-and-run). Run
`ctest --test-dir <build-dir> -R mcmc --output-on-failure`.

The test uses complex double tensors with virtual dimension 2 on four-site chains
and squares. An independent dense enumeration of virtual indices computes all
16 amplitudes. Graph contraction must match each within 1e-10 * (1 + abs(reference)) and leave consistent
site ownership metadata. Sampling then uses boundary bond cap 16, 30 BP steps,
1e-12 BP tolerance, and 1e-14 truncation thresholds. After 128 burn-in sweeps, all
16 state frequencies in 4096 samples must be within 0.04 absolute probability of
the dense reference, with total variation distance at most 0.10. A fixed seed makes each MPI configuration reproducible;
rank counts are not expected to produce identical sample sequences.

CTest runs with 1 and 2 ranks, one OpenMP/BLAS thread, and case-dependent 30/120/180-second timeouts.
The frequency check is a statistical smoke test, not a proof of unbiased sampling
on arbitrary graphs or a test of deliberate boundary truncation.
The CUDA build also provides these executables, without registering CPU CTest
launch commands for them. Explicit one-GPU and same-node two-GPU runs have passed
for chain/square distributions, manual/automatic heat-bath distributions, product
controls, fixed sites, and real/complex CLI/checkpoint cases. The heat-bath checks require
identical complete automatic/manual sample sequences within each rank count.
Other CPU regression modes are not all validated on CUDA. Large networks and
single/disconnected vertices remain outside this small-graph regression's scope.
See the separate 51-site probe for its limited GPU coverage.

## Diagonal and manual grouping tests

The manual CTest cases cover both diagonal orientations, an explicit singleton
group and three lines `{0}/{1,2}/{3}`, at 1 and 2 MPI ranks. Each case checks all
16 exact amplitudes and 4096 sample frequencies with the same thresholds above.
They use the workspace API. The CLI tests also cover
`--linefile`/`--groupfile` parsing and MPI distribution. CLI tests require
16 complex-valued-network samples to complete.

The test also rejects duplicate/missing/unknown members, an extra target-line
member, malformed/duplicate/missing input rows, and an invalid absorption path.
A direct singleton relabel check verifies ownership, IDs and edge-vector order.

`mcmc_regression diagonal03` and `diagonal12` select the two diagonal
orientations with automatic grouping. These cases use the original-graph
three-line planner and are included in the automatic grouping tests below.
Unsupported greedy partitions are rejected before boundary construction.

## Shared BMPS plan regression

For each of the 16 configurations in each numerical case, call both the tensor
and graph-only BMPS APIs directly and compare the resulting edge vectors. Complete
the remaining contraction and compare the scalar with independent dense enumeration
(tolerance 1e-10 * (1 + abs(reference))). This checks both automatic and explicit-group BMPS paths at
1/2 ranks. A diagonal case asserts identical nearest-unassigned automatic/manual schedules;
invalid/no-destination plans must throw before changing the supplied graph.

## Upper-region independence

The standalone eight-site upper/lower boundary test is part of the
[boundary-MPS suite](../boundary_mps/README.md).
Run the `bmps_upper_independence` executable or its CTest cases.

## Single-site end lines on a diagonal square

`mcmc_regression singleton-ends` uses edges `(0,1), (2,3), (0,2), (1,3),
(0,3)` and lines `{1}, {0,3}, {2}`. Both automatic grouping and manual groups
`{{{1,0,3}}, {{0,1,2},{3}}, {{2,0,3}}}` are tested, including a middle group
with neighbors on both sides and a singleton middle group. Each case compares
16 exact amplitudes and 4096 post-burn-in samples to the dense reference.
CTest runs this mode with MPI 1 and 2 ranks.

Environment initialization returns only the N-1 proper prefixes/suffixes of an
N-site line. All numerical cases directly check local output counts and site labels
in both directions; one-site calls must clear prefilled output buffers. The
singleton-ends mode also covers a three-site line `{0,1,2}` and a single-site line
`{3}` on a chain, so the first/interior/last update positions are exercised.

Boundary initialization checks also require the uppermost lower-boundary slot
(tensors, site IDs, messages and edge IDs) to be empty: no sampling update consumes
that final full-network boundary.

## Dimension-one sites

`mcmc_regression fixed-sites` tests a chain with lines `{0,1,2}/{3}` and
fixed-site masks 15 (all), 2 (middle site), and 5 (both ends of the three-site
line). Dense probabilities exclude states where a fixed site is 1. Mixed cases
use 4096 samples after 128 burn-in; the all-fixed case uses 64 after 4 burn-in.
CTest covers MPI 1/2 ranks and a CLI run with the default physical dimension 1.

All regression modes also check invalid line structures before numerical tests:
empty/single line, empty member line, duplicate/missing/unknown site, and a
nonadjacent-line edge. Both automatic/manual full APIs must throw invalid_argument
before tensor access; empty-batch workspace initialization must also reject and
preserve output containers.

## Absorption with long-bond relocation

`mcmc_regression shift-absorption` uses six sites and lines `{0,1}/{2,3,4}/{5}`.
Absorbing the lower line creates a 2--4 bond, which must move through site 3.
One graph has the existing 3--4 bond; another omits it and reverses the edge
vector and endpoint orientation. Explicit groups satisfy both sweep directions.
The second graph adds 3--5 so its middle-line environment suffix is connected.

All 64 physical states are checked against independent dense enumeration.
The direct tensor and graph-only increments must return identical ordered edges,
leave no long output-line bonds, and preserve output-site owners. Both left/right
environment initializations and both sampling sweep directions are exercised;
4096 samples are compared to the exact 64-state distribution. CTests run with
MPI 1/2 ranks and a 180-second timeout.

Standalone bond-shift tests additionally compare the exact ordered graph-only
result for real/complex, multihop, shared-neighbor and partial-input cases at
MPI 1/2/3 ranks; see [graph tests](../graph_tensor_contraction/README.md) and [BMPS tests](../boundary_mps/README.md).

## Automatic grouping

`mcmc_automatic_groups_1_rank` / `_2_rank` test grouping on the original graph.
Cases include a chain with lines `{0,1,2}/{3}` and a variant with an extra
0--2 edge. Automatic full-API sampling and
explicitly supplying its generated plan through the workspace API must produce
identical complete sample sequences with the same random seed. Both diagonal
orientations are also checked against their 16-state exact distributions.

Plan checks cover the 2x6 sparse-rung graph (rungs at positions 0, 2 and 4), exact
coverage and expected distance/tie choices, edge order/orientation invariance,
and nonmonotonic relabeling preserving the specified line order. Disconnected
neighbor components are rejected. A layout that this greedy policy cannot
partition is rejected even though a separately tested manual partition works.

The 12-site sparse-rung graph additionally runs through automatic workspace
initialization, exact numerical absorption, both environments and 64 samples
with physical dimension 1 and virtual dimension 2. This exercises the graph and
tensor paths but does not test a nontrivial 4096-state sampling distribution.
The initialized workspace must retain exactly the generated plan.

The policy and admissible-path distance are documented in
[automatic MCMC grouping](../../docs/automatic_mcmc_grouping.md).

## Heat-bath updates

`mcmc_regression heat-bath` checks method 1 on a correlated four-site state
against the enumerated exact distribution, and requires automatic/manual plans
to produce identical samples with the same seed. It also checks fixed sites,
a uniform binary product distribution (periodicity regression),
a three-state product distribution (including a zero-weight initial state),
zero/equal/extremely scaled amplitudes and rejection of unknown method IDs.
CTest runs it and the method-1 CLI smoke test with one and two MPI ranks.
Method 0 remains available explicitly; method 2 (full-state roulette) is now the default.

## Directory and targets

- `main.cc`: `mcmc_regression`, including grouping/line validation, environments,
  Metropolis/heat-bath statistics and end-to-end BMPS absorption integration.
- `data/`: graph, line and group fixtures used by the CLI tests.
- `CMakeLists.txt`: registration of this suite and the real/complex CLI targets
  built from `apps/mcmc_sampling`.
- [`heavyhex51/`](heavyhex51/README.md): 51-site graph, shallow circuit, manual and
  boundary-only automatic selection, independent product-state reference.

```sh
cmake --build build/tests --target mcmc_regression mcmc_sampling mcmc_sampling_complex
ctest --test-dir build/tests -R mcmc --output-on-failure
```

The 51-site experiments are optional executable runs, not part of default CTest.

## CLI workspace coverage

The CLI uses `MCMCWorkspace` for shared planning and per-batch caches. Manual
Metropolis and automatic-group heat-bath CLI cases run two batches at MPI 1/2
ranks and require the last sample of the second batch to be emitted.

## CLI input and checkpoint regression

`mcmc_options_regression` checks missing values for every value-taking option,
malformed/overflowing numbers, invalid ranges, comma-separated site IDs, help
precedence and valid boundary values.

`checkpoint_fixture.cc` writes a deterministic state `002` on unsorted edges
`(0,2), (0,1)`, with unequal virtual dimensions 2/3 and physical dimensions 1/1/3.
Real and complex CLI tests load it at MPI 1/2 ranks with two batches. At two ranks,
rank 0 owns no tensors. CTest fixtures generate the matching rank files before
sampling; the final sample must equal `002`.

CUDA builds provide the same real/complex fixture generators. Launch a generator
with one MPI process and one visible GPU; its `ranks` argument is the number of
checkpoint files to write for the subsequent sampler, not the generator's process
count. CPU and CUDA tensor serialization formats differ, so generate checkpoints
with the backend that will read them.

Fixed-site distribution regressions use pairwise heat bath. Near-balanced binary
Metropolis proposals can almost undo each other on the return sweep, so retained
samples can be strongly correlated even after burn-in. These cases keep the
same sample counts and distribution tolerances; Metropolis remains covered by
the default, manual, automatic-group, bond-shift and singleton regressions.


## Full-state roulette heat bath

`mcmc_full_heat_bath` enumerates all virtual indices independently for complex
D=2 K4 tensors (including reversed edge orientations), with physical dimensions
2/2/2/2 and 2/3/2/1. Every initial configuration, both directions and every
one-/two-site window are tested collectively using the distributed graph API.
All candidate absolute amplitudes and conditional probabilities must match
within 1e-11 relative/absolute tolerance. Every candidate prefix must stay on
the leading target owner and numerical caches must be rank-local. Move reuse
of actual candidate storage and element-wise exact immutability of all opposite-side
caches are checked. Candidate evaluation accepts no raw exterior network.

An exact roulette oracle uses a separate random stream for each owner and
checks repeated sweeps, overlapping pairs, both boundary sides, heterogeneous
dimensions and an empty rank zero. A fixed-charge 01/10 state must freeze under
one-site updates and visit both supported states under two-site updates. All
256 p16 pair candidates are checked without charge filtering, including when
rank zero owns no tensors. Scalar weight tests cover extreme scales, invalid
weights and roulette frequencies. Scalar-read failure on a nonzero owner is
reported once at the window end; a line-update failure aborts with a diagnostic.
Kernel tests run with CPU MPI 1/2/3 ranks.

`mcmc_regression full-heat-bath` checks chain/square/K4 distributions at both
window sizes using 128 burn-in round trips and 4096 retained samples per case.
Automatic and manual chain sequences must agree. It covers explicit-method and
default-method workspace APIs and the individual-container API. The existing
frequency thresholds (max error 0.04, TV 0.10) are statistical smoke tests.
Integration/CLI tests run on CPU MPI 1/2 ranks; CLI runs use ternary
states and two batches with the default method. Legacy Metropolis/Barker tests
remain registered. GPU, large-line performance and deliberately truncated-boundary
accuracy are not validated by this addition. Candidate contraction is restricted to window groups and cached sides; no
remaining-line snapshots or candidate-dependent exterior reconstruction are used.

An eight-group chain audit compares physical dimensions 1 and 3, window widths
1/2 and both directions. Opposite-side construction makes six cached graph
joins for one-site and five for two-site updates, once per sweep and independent
of candidate count; group preparation is called exactly eight times. Unused
pair suffix slot 1 must remain empty. Window joins are 14/42 for one-site
and 19/99 for two-site updates (d=1/3), each with exactly two cached region
inputs. Each selected leading intermediate is moved once. This checks cache
reuse work counts in addition to numerical correctness.


## Review-fix checks

Single shared legs are tested at chi=1/7/64 with no identity buffer, exactly one
binary graph contraction and two immutable-input clones at the cache-join entry.
A nonclosing join checks a shared leg following different external legs on both
operands, with all result components compared to an independent sum. Existing
K4 cases exercise multiple shared legs fused into a single composite bond. Source cache values remain
unchanged. An empty rank zero and distributed owners are included.

The eight-group production audit compares environment and candidate contraction
counts against the compiled native plans on each participating rank. It checks
zero identities, one environment construction, one selected-prefix move per
window, reused metadata, last-use prefix/suffix/right-candidate moves, and no line-wide shape reduction. Separate cache-helper
tests retain the two-input copy checks for the older standalone helper API.

Physical leaf IDs must have neither numerical tensors nor owners in the
prepared graph. Invalid initial states and mismatched virtual bond dimensions
must still terminate with a diagnostic after removing numerical physical leaves.


## Square-lattice contraction order

The kernel includes a three-column strip with separate target, upper and lower
rows, D=2 and chi=3. Independent enumeration of all 20736 original virtual-index
assignments gives the eight exact state amplitudes. Every candidate in every
one-/two-site window, all eight initial configurations and both directions are
compared with these sums (tolerance 1e-11 * (1 + reference)). Repeated production
sweeps are compared with independent owner-specific roulette streams. Both upper-
and lower-first member orders, split owners and an empty rank zero are covered.

An explicit trace checks the selected plan against its predicted MACs and
sizes. With two candidates the common boundary/side join runs once, followed
by target/other boundary per state: total MACs 612, largest output 36 elements.
With three candidates the common subtree is completed before target entry:
total MACs 756, largest output 72 elements. The two plans must give identical
values for identical raw states. Repeated metadata reuses the cached plan;
changing the candidate count requires another search.

The common intermediate stays on the owner of its smallest leaf ID (the side
in this fixture). Target entry and the final root retain the target owner.
Selected candidate storage is moved into the next proper environment. These
checks verify the constrained local plan, not wall-clock speed or a globally
optimal arbitrary-graph schedule.


## Honeycomb layouts

`mcmc_honeycomb_kernel_{1,2,3}_rank` checks a three-row/four-column brick-wall
patch with two hexagonal cycles, D=2 and boundary chi=3. Both alternating-rung
parities, reversed edges, indirect group members, all 16 target configurations,
widths 1/2 and both directions are covered. Independent sums enumerate all
93312 original virtual-index assignments. Owner-local candidates, immutable
opposite caches, storage moves, no identity buffers and exact roulette streams
are checked, including an empty rank zero. Traces exercise the selected common
and target-dependent stages. The planner may replace the old connectivity-first
order; no fixed entry order is required for all honeycomb orientations.

`mcmc_honeycomb_integration_{1,2}_rank` uses D=2 on the full patch with fixed
physical dimension-one boundary rows, leaving all 16 middle-row configurations.
Actual initialization and round trips, automatic/manual groups and both update
widths are tested with the existing 4096-sample thresholds (max error 0.04,
TV 0.10). Automatic/manual sample sequences must agree.

`mcmc_honeycomb_zigzag_{1,2,3}_rank` uses explicit unit-cell connectivity
A(m,n)--B(m,n), A(m,n)--B(m-1,n), A(m,n)--B(m,n-1), with three cells in each
of two rows. BFS from A(0..2,0) must return A0/B0/A1/B1 layers. It tests these
four automatic half-layers and two explicitly paired zigzag chains; the
production layer generator is unchanged and does not merge the pairs.
Six exterior physical sites have dimension one; every configuration of the
remaining six sites (64 states) is checked. Actual BMPS initialization with cap
16 is performed for every initial state, alternating automatic and explicit
grouping. All lines/windows, both widths and directions, 19456 candidate
amplitudes per MPI run and whole-line roulette results are compared with sums
on the original graph. BMPS normalization introduces a common exterior factor;
one scale per line is fixed once and required to match all candidates, windows,
widths and directions (tolerance 1e-10). Storage moves, split owners and an empty
rank zero are included. Adjacent positions in a separate A/B layer need not be
joined by an original bond; adjacent sites in the paired zigzag chains are.

These are finite-patch CPU checks, not GPU, large-scale performance or general
finite-chi accuracy tests. A trivalent graph obtained by tensor factorization
need not have honeycomb topology; ancilla transformations are not implemented
or verified by these fixtures.

`mcmc_regression zigzag-distribution automatic|paired 1|2` adds distribution
checks on the same 12-site unit-cell patch, with six fixed physical dimension-one
sites and six binary sites (64 allowed configurations). Both automatic A/B
layers and paired zigzag lines are tested at both update widths on CPU MPI 1/2
ranks. Each case discards 128 round trips and compares 4096 retained samples
with normalized Born probabilities from independent enumeration of the original
virtual indices. The existing thresholds remain max absolute error 0.04 and
total variation 0.10. Boundary-MPS bond cap is 16. These are statistical
regressions, not mixing-time measurements or new GPU distribution checks.


## Move-only cache and commit fixes

Numerical boundary caches are reset with clear/resize, without requiring Tensor
copy construction (CUDA Tensor is move-only). Commit takes prepared raw target
slices and moves the selected slice into its original owner-local site slot;
it accepts neither a tensor context nor W and performs no physical projection.

The actual commit template is instantiated with std::unique_ptr<int> in every
kernel mode. Pointer identity, consumed selected slices, preserved unselected
slices, off-line states, permuted local site/target orders and an empty rank
zero are checked. This also rejects a regression to copying move-only values.
Existing numerical tests check original-axis tensor values after the same
production commit for both widths and directions, including mixed physical
dimensions, scalar singleton tensors and overlapping pairs. Separate cached-join
checks cover unequal virtual-leg dimensions and axis order.


With no side environment, the first boundary member is borrowed directly from
its immutable owner-local cache. Only an actual pre-target contraction creates
an owned common intermediate. The graph join still clones each immutable input
once internally. The source boundary remains alive through all candidate states;
no separate numerical seed copy is needed. This preserves absorption order,
owner placement, opposite environments and raw-slice commit.

The heat-bath kernel also checks planned group contraction: common-subtree reuse
for two/three candidates, line-local plan reuse, target-owner retention,
prepared multi-bond fusion against a direct complex sum, and fatal failure
cases in separate bounded MPI jobs (including the no-step singleton case).
PMPI counts require no Allreduce in group/candidate/graph execution and one
scalar-failure notification per window. The production line has no Allreduce. Direct line calls share candidate counts
only on each window's communicator; prepared sampling calls already know them.
Cross-owner contraction dimensions are checked by the native graph operation. Ten fatal cases
per MPI configuration also cover invalid initial state, cross-owner dimension
mismatch, missing raw/commit tensors and an all-zero window.
Existing square, brick-wall and zigzag exact-candidate tests cover both window
widths and sweep directions with distributed owners.

The production line entry point is also checked with parent ranks 0 and 2 as
participants: rank 0 owns boundaries only and rank 2 owns targets. Rank 1 waits
outside the function until rank 2 completes, proving no parent-wide collective
is needed, including communicator creation. Both widths and directions are
compared with the exact sweep oracle. PMPI counts require zero parent
Allreduce/Bcast calls, zero line shape reductions, and no
communication for nonparticipants. Original graph/owner metadata, unrelated
tensors/states and nonparticipant RNG are preserved.

Sampling-call communicator reuse is checked separately from contraction-plan
selection. The three-rank kernel reuses one {0,2} communicator across widths and
directions while rank 1 waits outside. PMPI requires one creation/free pair.
The sampling regression also counts one creation/free per unique participating
owner set across every sample and both sweeps. Direct line calls without an
external registry retain their temporary-communicator behavior.


## Original-edge heat-bath execution

Default production sweeps use `mcmc_edge_plan.h` and `mcmc_edge_heat_bath.h`.
Square and brick-wall production oracles exercise baseline and both D < chi
and D > chi model choices, one/two sites, both directions and distributed
owners. Zigzag production oracles use model-selected native plans for automatic
A/B and paired layers. Earlier standalone cached-helper checks remain as
independent compatibility/reference tests; their small-graph work counts do
not describe the new native execution path.

A native K4-region fixture uses interleaved edges and unequal factors. A common
subtree changes the target and an exterior neighbor's axes. All three candidate
output tensors are compared componentwise with independent sums; the target
owner and immutable raw caches are checked. The common contraction runs once,
the dependent contraction three times, with ten total private input copies
across owners and no intermediate-step recopy. The three-rank case includes an
owner holding only the external neighbor. PMPI checks zero executor collectives.

The 4096-sample integrations continue to compare exact Born probabilities for
square/brick-wall and automatic/paired zigzag layouts, both window widths.
Communicator creation/free counts cover the full sampling call, with the exact
native-plan owner sets rather than a single three-line owner set. Original
physical dimensions/model D are shared once at sampling entry; no shape
collective or planner broadcast occurs inside a line.

A four-row metadata fixture checks the boundary checkpoint model against the
existing edges-only BMPS replay: unabsorbed outer-row bonds remain D and only
absorbed/truncated boundary bonds become chi.

The native scalar-join regression prepares one descriptor for six state pairs on
a three-factor cut (2 x 3 x 4), compares complex scalar values with an independent
sum, checks immutable A candidates and last-use B moves, and verifies that no
collective is introduced. It runs with 1/2/3 ranks, including an idle rank.

## Alternating two-site layouts

`mcmc_alternating` checks the batched wrapper against explicit initialization and
sampling for `{1,1}` and `{2,3}` round trips, including final RNG equality and zero
outputs. A six-site rectangle checks switching between two and three lines
with automatic/manual groups and nontrivial boundary construction. A four-site square with exactly one particle tests conservation of total
charge, row-charge confinement under horizontal-only updates, and transfer between
rows under alternating updates. After 128 burn-in cycles, 4096 samples are compared
with the analytic Born probabilities (maximum error < 0.04, TV < 0.10).
CTest runs on 1/2/3 ranks, with an empty third rank; complex CLI tests run on 1/2
ranks with manual first-layout and automatic second-layout grouping. Option tests
cover incompatible modes and invalid phase counts; fatal-input checks cover
invalid second layouts/groups, including zero-output requests. These tests do not establish
accuracy at arbitrary boundary truncation or on arbitrary layouts.

## Initial-state files

`mcmc_initial_states` verifies arbitrary file/local site order, batch mapping,
empty owners and mixed physical dimensions (including dimension-one ancillas).
`mcmc_initial_states_cli_*` runs on 1/2/3 ranks, including empty rank zero, with
16 physical states and two batches. Singleton lines with two-site updates keep
states unchanged so both single and alternating CLI output can be checked exactly.
Malformed/missing files, missing/duplicate/unknown sites, wrong batch columns,
negative/noninteger/overflowing/out-of-range states must terminate with diagnostics,
even with zero requested samples. Charge and nonzero Born support are not inferred.
