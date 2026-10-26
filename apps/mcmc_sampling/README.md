# Markov-chain Monte Carlo sampler

This application samples a tensor-network state using general-graph boundary MPS
and the TCAPI tensor interface.

## Portable checkpoint input and output

See [portable TPS checkpoints](../../docs/portable_checkpoints.md) for independent
site/message files shared by estimator and sampler, including restart with a
new MPI process count. This is a separate format from the legacy rank-local
options described below. No external library dependency is added.

## Build

Configure the [common test entry point](../../test/README.md#build-and-run)
for the CPU or CUDA backend, then build. The CPU configuration needs MPI, gqten,
HPTT, OpenMP and BLAS/LAPACK; see [CUDA build selection](#cuda-build-selection)
for the GPU configuration:

```sh
cmake --build <build-dir> --target mcmc_sampling mcmc_sampling_complex
```

`mcmc_sampling` uses double; `mcmc_sampling_complex` uses complex double.
Use the complex executable for complex checkpoints from estimator applications.
Checkpoint backend, scalar type and MPI partition must match the producing run.

## Example

Create `chain.txt` containing `0 1`, `1 2`, `2 3` on separate lines, then run:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpirun -np 1 <build-dir>/mcmc_sampling \
  --inputedge chain.txt --init_line 0 --p_bond_dim 2 --v_bond_dim 1 \
  --num_samples 100 --num_batches 1
```

`--inputedge` supplies undirected edges. `--init_line` supplies comma-separated
initial boundary sites; `--linefile` instead supplies one line of site IDs per
boundary. `--inputname` loads rank-specific TPS files in SaveTPS format; otherwise
tensors are random. Results are printed as
`sample <sample>.<batch>: <bitstring>`, in ascending site-label order.

## Options

Run `<build-dir>/mcmc_sampling --help` (or `-h`) to display all options without
loading inputs. Under MPI, only rank 0 prints help and all ranks exit successfully.
The real and complex executables accept the same options.

| Option | Default | Description |
|---|---|---|
| `-h, --help` | `—` | Show this help and exit successfully; no input files are required. |
| `--inputedge FILE` | `required` | Undirected graph edge file, one pair of site IDs per row. |
| `--portable_input DIRECTORY` | `none` | Load shared portable TPS and its graph on the current MPI partition. |
| `--initial_statefile FILE` | `none` | Initial physical basis indices, one column per batch; omitted means random. |
| `--inputname PREFIX` | `none` | Load rank-specific TPS files; otherwise generate random tensors. |
| `--outputname NAME` | `none` | Displayed only; does not write a file. Redirect stdout to save samples. |
| `--seed N` | `1729` | Random seed for tensor generation, initial states and sampling. |
| `--p_bond_dim N` | `1` | Physical dimension of random tensors; 1 gives fixed state zero. |
| `--v_bond_dim N` | `1` | Virtual bond dimension of random tensors. |
| `--init_line IDS` | `none` | Comma-separated initial boundary IDs for automatic line construction. |
| `--linefile FILE` | `none` | Explicit lines: one row of space-separated site IDs per line. Overrides --init_line. |
| `--groupfile FILE` | `none` | Manual three-line groups; omitted means automatic grouping on the original graph. |
| `--line_mode single\|alternating` | `single` | Select one layout or alternate two layouts (method 2, update_sites 2). |
| `--second_linefile FILE` | `none` | Required second layout in alternating mode. |
| `--second_groupfile FILE` | `none` | Optional second layout groups; otherwise automatic. |
| `--first_round_trips N` | `1` | Positive number of first-layout round trips per alternating cycle. |
| `--second_round_trips N` | `1` | Positive number of second-layout round trips per alternating cycle. |
| `--bmps_bond_dim N` | `1` | Maximum retained boundary-MPS bond dimension. |
| `--bmps_bp_tolerance X` | `1e-8` | Belief-propagation residual stopping tolerance. |
| `--bmps_max_bp_steps N` | `10` | Maximum BP iterations per boundary truncation. |
| `--bmps_tg_err X` | `1e-8` | Target boundary truncation error. |
| `--bmps_sv_min X` | `1e-8` | Singular-value cutoff for boundary truncation. |
| `--num_samples N` | `1` | Saved samples per batch; one round trip in single mode, both phases in alternating mode. |
| `--num_batches N` | `1` | Number of chains initialized and sampled. |
| `--sampling_method N` | `2` | 0: Metropolis; 1: pairwise Barker; 2: full-state roulette heat bath. |
| `--update_sites N` | `1` | One or two consecutive line sites; two requires method 2. |

Supply `--inputedge` and either `--init_line` or `--linefile` for sampling.
If both line options are present, `--linefile` takes precedence. The physical and
virtual random-tensor dimensions do not apply when loading `--inputname`.
Dimensions and batch count must be positive integers; sample count and maximum
BP steps may be zero. These values are limited to the signed `int` range. The
seed is an unsigned 32-bit integer. Tolerances must be finite and nonnegative.
Missing values, malformed numbers, out-of-range values and unknown options are
rejected with an error. `--help` takes precedence over other arguments.

## Sampling methods and output

`--sampling_method 2` (default) uses full-state roulette heat-bath updates.
The default `--update_sites 1` considers every state of one site. Set
`--update_sites 2` to consider all state pairs of consecutive sites in the line,
without charge filtering. The current state is included, and candidates are
selected with probability proportional to their Born weight.

Two-site windows overlap: `(0,1), (1,2), ...`, with the order reversed on the
return pass. Singleton lines are skipped in two-site mode. To enable it, add:

```sh
--update_sites 2
```

Methods `0` (Metropolis) and `1` (two-state Barker) remain available explicitly
and require `--update_sites 1`. Method 1 is full heat bath for binary physical
dimensions. Full-state heat bath allows staying in place and avoids the
always-flipping Metropolis periodicity for uniform binary states.

Boundary-MPS truncation remains approximate. All-zero or nonfinite candidate
weights are errors. Two-site updates alone do not guarantee access to the full
support; in particular, fixed-total-charge states cannot transfer charge between
lines with a fixed line layout. `--update_sites` does not change that layout.

Samples are written to standard output; redirect it to a file to save them.
`--outputname` does not create a file. In single mode each saved sample follows one full
forward/backward sweep; alternating mode saves after both phases finish. Burn-in removal and thinning must be done by the caller.
With the default physical dimension of one, all reported states are zero.
For physical dimensions above two, the output contains state indices rather than
binary digits; indices are concatenated without separators.

## Alternating line layouts

To exchange along two directions, supply two valid layouts of the same graph:

```sh
--line_mode alternating --update_sites 2 \
--linefile horizontal.lines --second_linefile vertical.lines \
--first_round_trips 2 --second_round_trips 3
```

One output sample is the state after two full first-layout round trips **and**
three full second-layout round trips in this example. Intermediate states are
not output. Both counts default to one. Swap the files to reverse the phase order.
Each phase continues from the latest state with the same rank-local RNG and
rebuilds its boundary MPS. Plans are reused across batches and round trips within
that phase. `--second_groupfile` optionally specifies the second layout's groups.
Phase-specific options require alternating mode; method 2 and update_sites 2 are
required. Zero output samples validate the layouts without running the phases.

Alternating layouts can permit charge transfer between rows, but do not guarantee
access to every state. For fixed-charge tensors, the initial configuration must
have nonzero weight; use `--initial_statefile` to supply an appropriate configuration.

## Initial physical configurations

Use `--initial_statefile states.txt --num_batches 2`, for example, with:

```text
# site_id  batch0  batch1
0          3       12
1          5        6
2          0        0
```

Each row gives a site ID followed by exactly one physical basis index per batch.
Row order is arbitrary; blank lines and `#` comments are allowed. Every graph site,
including ancillas, must occur exactly once. Indices must satisfy `0 <= state < d_i`
for that site's physical dimension (a dimension-one site must be zero).
Missing/duplicate/unknown sites, malformed integers, incorrect column counts and
out-of-range states are errors, including when `--num_samples 0`.

This input works with either line mode and with random or checkpoint tensors.
Indices refer to the tensor's physical basis order; the sampler does not interpret
blocked-site charges or check that the initial configuration has nonzero Born weight.
The file is read and shared once at startup, then each owner keeps its local states.
No random draws are used for file-based physical initialization. Without this option,
the previous random initialization is retained.

## Checkpoint input and edge order

When using `--inputname PREFIX`, provide `--inputedge FILE` in the **same edge
order used when saving the tensors**. Virtual tensor axes follow incident-edge
order. The sampler preserves the file order and rejects duplicate undirected
edges instead of sorting them. The TPS format does not store the graph, so an
incorrect edge order cannot in general be detected when dimensions coincide.

Each rank loads its own TPS file, including valid files containing no local
tensors. Without `--initial_statefile`, initial physical states are drawn using each local
tensor's last-axis dimension; physical dimensions may differ between sites. Empty ranks keep empty
state arrays and still participate in the collective sampling calls. Backend,
scalar type and MPI ownership must match the saved data. A missing checkpoint
file is reported before attempting to deserialize it.

## Lines and environment groups

Use `--init_line` to construct lines from a starting boundary, or `--linefile`
to specify every line explicitly. Within each line, site order determines the
left-to-right sampling sweep.

Without `--groupfile`, the sampler constructs groups automatically from the
original graph. For each target site it first reserves an available directly
connected neighboring-line site, then assigns the remaining sites by growing
connected groups in distance order. See the
[automatic grouping policy](../../docs/automatic_mcmc_grouping.md) for tie-breaking
and path restrictions. If automatic grouping rejects a layout, specify manual
groups or change the lines.

Automatic and manual groups are fixed at initialization and reused for upper/lower
boundary absorption and left/right environment construction. Groups are not
selected again during sampling.

### Manual groups

Use `--linefile lines.txt --groupfile groups.txt` to supply explicit lines and
groups. Each group is anchored at one target site and contains its assigned sites
from the immediately neighboring lines.

For lines `0 1` and `2 3`, a paired grouping file is:

```text
# zero-based target-line index, target site, neighboring-line sites
0 0 2
0 1 3
1 2 0
1 3 1
```

The first row describes group `{0,2}` for target site 0 on line 0. The last two
rows explicitly define absorption in the opposite direction. For a middle line,
include the assigned sites from both neighboring lines. Each target site needs
exactly one row; the target is implicitly a member, so do not repeat it. A row
with only line index and target site specifies a singleton group.

Example with a singleton:

```text
0 0 2 3
0 1
1 2 0
1 3 1
```

Within each target line, all neighboring-line sites must be assigned exactly once.
The neighbor order on each row is the boundary absorption order (filtering to the
side being absorbed). This order must permit successive contractions into the
target; no automatic reassignment is performed. Group rows may be in any order;
group order follows target-site order in `--linefile`.

An example and input files are in `test/mcmc_sampling/data`:

```sh
mpirun -np 2 <build-dir>/mcmc_sampling_complex \
  --inputedge test/mcmc_sampling/data/diagonal03.edges \
  --linefile test/mcmc_sampling/data/two.lines \
  --groupfile test/mcmc_sampling/data/paired.groups \
  --p_bond_dim 2 --v_bond_dim 2 --bmps_bond_dim 16 --num_samples 16
```

### Supported layouts

Both automatic and manual grouping require:

- At least two nonempty lines covering every graph site exactly once.
- Edges within a line or between adjacent lines; no self-loops or duplicate edges.
- Every neighboring-line site assigned exactly once for each target line.
- Connected groups and connected growing prefixes/suffixes of groups. Singleton
  groups are allowed.
- A member order that permits successive boundary absorption into the target.

The sampler checks the grouping before constructing boundary tensors. Automatic
grouping is a greedy procedure and may reject a layout even when a suitable
manual partition exists. The root rank reads manual groups and distributes them
to all ranks. Input errors terminate the MPI job with an error message.

## Library interface

Set `MCMCWorkspace::lines` and, optionally, `manual_grouping`, then call
`mcmc_initialize` followed by `mcmc_sampling`. Leave `manual_grouping` empty to
select automatic grouping. Keep the graph, line order and grouping consistent
between initialization and sampling; reinitialize after changing them.

`mcmc_sampling_alternating` accepts initial batched states, two layouts/groupings
and two positive round-trip counts. The output vector's length sets the number
of complete cycles; its entries are filled by the function. Empty grouping entries
select automatic grouping. No caller-owned workspace is required.

For the individual-container overloads, pass a `tnbp::MCMCGrouping<IntT>` pointer
as the final optional argument to initialization and sampling. Its structure is
`groups[line][target position][member]`, with the target as the first member.
Physical dimensions may differ between sites, including dimension-one fixed sites.
See [graph and tensor conventions](../../docs/graph_tensor_conventions.md) for
edge order and ownership requirements.

## Validation and accuracy

[Sampling regression tests](../../test/mcmc_sampling/README.md) cover small graphs,
automatic/manual groups, both update methods, fixed sites, CLI input validation
and real/complex checkpoints including empty MPI ranks. A
[51-site example](../../test/mcmc_sampling/heavyhex51/README.md) provides a larger
graph with manual and automatic grouping.

Boundary-MPS truncation can make the sampled distribution approximate. These
checks do not guarantee accuracy for arbitrary graphs or parameter choices;
check convergence with boundary bond dimension and sampling length for your problem.

## CUDA build selection

The common test build accepts `-DTNBP_TEST_BACKEND=cuda` and builds both
`mcmc_sampling` (real double) and `mcmc_sampling_complex` (complex double).
See [backend configuration](../../test/backend_integration/README.md#build)
for CUDA and cuTENSOR paths. The CLI selects `tcapi::cuda::Tensor` with
`USE_TCAPI_CUDA`; CPU remains the default. Each MPI process must see exactly
one GPU, selected before any TCAPI context is created.

Real/complex CLI runs, native CUDA checkpoint loading, small exact-distribution
comparisons, and 51-site manual/automatic grouping have been exercised on one
GPU and two GPUs on the same node. The 51-site GPU runs use 64 retained heat-bath
samples: product-state samples match an independent scalar reference; the D=2
cases check execution and finite output, not distributional accuracy.
CPU and CUDA tensor checkpoint formats are not interchangeable. Load checkpoints
written by the same backend, with matching scalar type and MPI ownership.

## Site ownership

Use `--ownerfile FILE` for `site_id mpi_rank` rows; omitted selects the
automatic contiguous partition. The selected assignment is printed once as
`# site ownership: ...` and `# MPI rank N: ...` without timestamps. See
[site ownership](../../docs/site_ownership.md) for input validation, direct
initialization/loading and legacy checkpoint restrictions.
