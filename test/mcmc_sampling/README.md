# Sampling migration regression

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
on arbitrary graphs or a test of deliberate boundary truncation. GPU sampling,
large networks and single/disconnected vertices remain outside
this initial test scope.

## Diagonal and manual grouping tests

The manual CTest cases cover both diagonal orientations, an explicit singleton
group and three lines `{0}/{1,2}/{3}`, at 1 and 2 MPI ranks. Each case checks all
16 exact amplitudes and 4096 sample frequencies with the same thresholds above.
They use the workspace API, while the CLI tests exercise the full API and the
`--linefile`/`--groupfile` parsing and MPI distribution paths. CLI tests require
16 complex-valued-network samples to complete.

The test also rejects duplicate/missing/unknown members, an extra target-line
member, malformed/duplicate/missing input rows, and an invalid absorption path.
A direct singleton relabel check verifies ownership, IDs and edge-vector order.

`mcmc_regression diagonal03` and `diagonal12` select the earlier automatic
reproducers. Their original failures are archived in the workspace experiment.
These automatic cases now use the original-graph three-line planner and are
included in the automatic grouping regression below. Unsupported greedy
partitions are rejected before boundary construction.

## Shared BMPS plan regression

For each of the 16 configurations in each numerical case, call both the tensor
and graph-only BMPS APIs directly and compare the resulting edge vectors. Complete
the remaining contraction and compare the scalar with independent dense enumeration
(tolerance 1e-10 * (1 + abs(reference))). This checks both automatic and explicit-group BMPS paths at
1/2 ranks. A diagonal case asserts identical nearest-unassigned automatic/manual schedules;
invalid/no-destination plans must throw before changing the supplied graph.

## Upper-region independence

The standalone eight-site upper/lower boundary test has moved to
[`boundary_mps/upper_independence.cc`](../boundary_mps/README.md).
Run the `bmps_upper_independence` executable or its CTest cases;
`mcmc_regression upper-independence` is no longer an entry point.

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

## Fixed automatic grouping

`mcmc_automatic_groups_1_rank` / `_2_rank` test the new original-graph policy.
The chain `{0,1,2}/{3}` and its extra 0--2 edge reproduce the previous unsafe
indexing and premature target consumption. Automatic full-API sampling and
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
Method 0 remains the default Metropolis kernel.

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
