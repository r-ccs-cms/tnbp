# Site ownership in applications

`apps/estimator`, `apps/estimator_restart`, and `apps/mcmc_sampling`
(`mcmc_sampling` and `mcmc_sampling_complex`) accept `--ownerfile FILE`.
Each row assigns one site to a rank in the run's communicator:

```text
# site_id mpi_rank
0 0
1 1
2 1
3 0
```

Rows may be unordered. Blank lines and `#` comments, including trailing comments,
are allowed. Exactly two decimal integers are required per data row. Every graph
site must occur once, with no extra sites; ranks must be in `[0, MPI_size)`.
Rank zero reads the file, validates it and broadcasts the selected owners. A
bad file is reported to all ranks before tensors are initialized or loaded.
Portable input still requires every rank to own at least one site. Fresh
initialization may have empty ranks.

Without `--ownerfile`, the apps choose a balanced contiguous partition of sorted
site IDs using `get_range`: lower ranks receive one extra site when division
has a remainder. Fresh initialization and portable input use the same rule.
This unifies portable's former `floor(site_position * MPI_size / site_count)`
rule with the existing product/random-state initialization convention.

The actual assignment is printed once, by rank zero, without a timestamp:

```text
# site ownership: owners.txt
# MPI rank 0: 0 3
# MPI rank 1: 1 2
```

The source is `automatic` for the default partition, the file path for an
explicit assignment, or `checkpoint` for legacy input without an ownerfile.
Empty ranks have an empty site list. The output records the final assignment,
not just the requested option.

## Initialization and portable loading

The app chooses ownership before tensor creation. Product-state and random-state
initializers construct tensors directly on their assigned ranks; no subsequent
redistribution is performed. For portable input the same map is passed to
`LoadTPSPortable`, so each rank directly reads its assigned V files and, for
estimator restart, the required E files. Sampling still reads V only.

For example, load a prepared state using a new two-rank assignment:

```sh
mpirun -np 2 ./estimator --circuit evolution.qasm --sparse_pauli pauli.txt \
  --portable_load prepared-state --ownerfile owners.txt --step_end 10

mpirun -np 2 ./mcmc_sampling_complex --portable_input prepared-state \
  --ownerfile owners.txt --linefile lines.txt --num_samples 1000
```

The estimator executable above is built from `apps/estimator_restart`.
Ownership is a run-time choice and is not stored as part of the portable state.
The writer MPI size and placement do not constrain readers: for example, a
checkpoint saved with two ranks can be read with four ranks and a newly chosen
ownerfile using ranks 0 through 3. Only the current communicator is checked.

## Legacy rank-local input

Legacy `--loadname` / `--inputname` files already contain the ownership map and
local tensors. Without `--ownerfile`, that saved map is used and printed. A
supplied ownerfile must match it exactly; changing the assignment requires
portable input. The apps reject a mismatch rather than ignoring the option.

## Initialization API

`InitTensorProductState` now takes its existing `Site_To_MpiRank` parameter as a
**const input map**. It must be complete, valid, and identical on all ranks.
The initializer no longer assigns owners or reconstructs the map with an
Allreduce. Choose a default explicitly before calling:

```cpp
auto sites = tnbp::GetSiteIndexFromBond(edges);
auto owners = tnbp::ContiguousSiteOwners(sites, comm);
tnbp::InitTensorProductState(ctx, edges, pdim, V, SiteIdx, owners, E, EdgeIdx, comm);
```

`tnbp::Init(ctx, pdim, edges, owners, comm, state, E, EdgeIdx)` likewise requires
an explicit map. `ValidateSiteOwners` checks coverage, rank ranges and agreement
between ranks. Tensor axis order remains defined by the ordered graph, not the
owner map. The initializers keep no fallback partition implementation.
Ownership-file parsing and comment-prefixed logging live in the shared app
helper `apps/site_ownership.h`, not the tensor library.

## GPU validation

ROQUO GB200 validation passed on one GPU shared by 1/2/4 MPI ranks and on two
GPUs with 1/2 ranks. All 50 application ownership acceptance cases passed,
including fresh product/random initialization, custom direct reads, empty ranks
for fresh states, invalid-file rejection and comment-prefixed assignment logs.
A checkpoint written with two ranks was read with a custom four-rank map and
preserved all V/E payload bytes. A later job on another node also loaded that
state with a new owner map on two GPUs and preserved V/E bytes. Sampling was
checked against the known cluster-state distribution. This does not validate
simultaneous multi-node execution or large-state performance.

The completed checkpoint format and string/POSIX path implementation were
subsequently validated on ROQUO (92 cases). Custom four-rank reads of a
checkpoint saved by two ranks preserved all V/E bytes; a later two-GPU job on
another node did so with a different, uneven owner map. Sparse site IDs
`{0,1,9,1000000}` and signed IDs `{-4,2,7,20}` also passed roundtrips and
custom-owner sampling. Ownership file keys are actual labels, not dense site
positions; the supplied line file must use the same labels.
