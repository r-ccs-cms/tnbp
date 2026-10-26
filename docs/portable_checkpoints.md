# Portable TPS checkpoints

`estimator_restart` can save a TPS that `mcmc_sampling_complex` reads directly
with a different MPI process count. No converter, root Tensor gather, or Tensor
redistribution is used. Estimator restart loads its owned vertices and both directed messages on incident
edges. Sampling loads only its owned vertices and never opens, hashes or decodes
E payload files. Both use TCAPI file I/O on the shared filesystem.

## Save with estimator and sample on a new partition

Build `apps/estimator_restart` with the usual backend configuration, or build
the `estimator_portable` CMake test target (the same application source).
Save after complete circuit cycles:

```sh
mpirun -np 2 ./estimator --circuit circuit.qasm --sparse_pauli observables.txt \
  --step_end 10 --portable_save /shared/run/checkpoint-10
```

The parent directory must exist. The checkpoint directory must be new. Then,
using the matching backend and scalar type:

```sh
mpirun -np 4 ./mcmc_sampling_complex \
  --portable_input /shared/run/checkpoint-10 --linefile lines.txt \
  --bmps_bond_dim 16 --sampling_method 1 --num_samples 1000
```

The graph is read from the checkpoint; `--inputedge` is optional. If supplied,
it must match the saved graph and edge order. Lines/groups remain sampler
settings and are constructed against that graph. Real checkpoints use the real
sampler, complex checkpoints use the complex sampler. `--portable_input` and
legacy `--inputname` are mutually exclusive.

The library checkpoint contains only the tensor state, ordered graph and
format/integrity information. Sampling prints the source state digest and never
reads estimator restart information.
Portable sampling seeds independent rank streams from `(seed, rank)`; changing
MPI count therefore changes the random sequence, while retaining the target
state. Legacy sampling retains its previous seed behavior.

To continue estimator evolution instead:

```sh
mpirun -np 4 ./estimator --circuit circuit.qasm --sparse_pauli observables.txt \
  --portable_load /shared/run/checkpoint-10 --step_end 20 \
  --portable_save /shared/run/checkpoint-20
```

Estimator writes its own optional `estimator.resume` file alongside the state
files. It contains the circuit hash, next step, BP/truncation settings and
complete-cycle phase, bound to the state digest and protected by its own checksum.
Only the estimator app reads/writes this file and decides how to resume.
If absent, the state is an initial state and estimator starts at step zero.
If present but corrupt, unsupported or bound to another state, estimator rejects
it; sampler ignores it. A leftover `estimator.resume.pending` is rejected by
estimator, rather than treated as ordinary state-only input.

Estimator owns an outer save transaction. It reserves
`<destination>.estimator.pending/`, asks the unchanged state library to save into
its `state/` child, and writes its resume file there. Only after both succeed
does rank zero rename that complete child directory to the final destination.
Thus readers of the final name see the state and resume information together.
A failure leaves the staging directory for inspection and does not publish the
final name. Existing destinations are rejected. The generic state library
continues to save only tensors/graph and never manages estimator metadata.
This is atomic publication via same-filesystem directory rename, not fsync
power-loss durability.

With a valid resume file, estimator validates its version and phase. Circuit
identity is the SHA256 of the exact QASM file contents (including comments and
whitespace). With the same circuit, BP/truncation settings must match and the
saved next step is inferred. With a different circuit, the saved V/E state is
used as the initial state of the new circuit, starting at step zero; new
BP/truncation settings may be supplied. An explicit `--step_start` must equal
that inferred start (saved next step or zero).

The saved graph and its edge/axis order are used to construct the new circuit's
TPO, including for a circuit containing only single-site gates. Circuit sites
must belong to the saved graph; the graph is not rebuilt from the new circuit.
An explicit non-default topology backend must match the saved graph.
Disconnected graphs are valid, but each multi-site gate must have paths between
its operands in that graph. An unrepresentable gate is rejected by TPO
construction instead of being silently omitted. Gates within a component and
single-site gates remain valid; graph edges/axis order are never extended.

For example, apply a Trotter circuit ten times after state preparation:

```sh
mpirun -np 4 ./estimator --circuit trotter.qasm --sparse_pauli observables.txt \
  --portable_load /shared/run/prepared-state --step_end 10 \
  --portable_save /shared/run/trotter-10
```

When `trotter.qasm` differs from the preparation circuit, this executes steps
`[0,10)`. Loading `trotter-10` with the same Trotter QASM and `--step_end 20`
then executes `[10,20)`. The estimator resume file records the new circuit and its
next step. Set `--step_end` to the exclusive endpoint. Loading the same circuit
with `step_end == next_step` saves without evolution.
Legacy rank-local `--loadname`/`--savename` continue to use their existing format;
legacy-to-portable conversion is not implemented here.

## Layout and integrity

```text
checkpoint-10/
  manifest.txt
  COMMITTED
  vertices/site-000000.ten
  vertices/site-000001.ten
  ...
  messages/edge-000000-forward.ten
  messages/edge-000000-reverse.ten
  ...
  estimator.resume             # optional; owned by estimator app
```

`manifest.txt` starts with `TNBP_PORTABLE_TPS 3`. It contains the codec/ABI,
axis convention, sorted site IDs, ordered edges with stable snapshot edge IDs,
and file records. File records contain path, length,
SHA256 and shape. No circuit, step, estimator settings, algorithm/residual or
application key/value fields are stored in this format. All counts are checked
decimal integers; strings use C++ `std::quoted`.
No JSON, external parser, hash command, or new library dependency is used.
SHA256 is implemented incrementally in the library and checked against known
vectors by the regression suite.

Vertex virtual axes follow the saved incident-edge order; the single physical
axis is last. For saved edge `(a,b)`, `forward` means `a` to `b`, `reverse` means
`b` to `a`. E axis zero is ket, axis one is bra. Directions are independent
records. The loader never sorts/reverses graph edges or transforms tensors.
The state digest covers the state format, codec, axis convention, graph and
logical record hashes; it excludes MPI ownership and the estimator
resume file.

Writers use a sibling `<destination>.pending` directory. Each V has one owner
writer, each directed E has its source-site owner writer. Both boundary copies
must serialize identically; nonwriters use temporary check files in pending,
which are removed before publication. All ranks agree on write success and
record coverage before the manifest and its SHA256 `COMMITTED` marker are
written and pending is renamed to the final directory. An incomplete pending
directory is never an input checkpoint. It is retained on failure for inspection;
choose a fresh destination for retry. Existing checkpoints are not overwritten.
This provides publication at directory granularity on filesystems with atomic
same-filesystem rename, not power-loss durability via fsync. Each destination is
one immutable generation; no automatic `latest` pointer is maintained.

Readers check manifest integrity, owner coverage, graph/shape/codec consistency,
and only the selected payload files' size, SHA256 and native headers before
allocation. The default per-record allocation limit is 64 GiB and is adjustable
in the library API. Errors are agreed between live ranks before returning;
process failure still requires MPI/launcher abort or timeout. Loading never
normalizes, fixes phases, runs BP, truncates, or applies state preparation.

## Library API and limits

`ReadTPSCheckpoint(directory)` returns the graph and logical tensor records.
`PortableContiguousOwners(sites,comm)` supplies the CLI's
balanced contiguous partition (lower ranks receive remainder sites).
The applications accept `--ownerfile FILE` to override it and print the actual
assignment; see [site ownership](site_ownership.md). `LoadTPSPortable<T>(ctx,directory,owners,comm)`
also accepts an explicit owner map and returns local V/E, SiteIdx/EdgeIdx,
new owners, communication order and descriptor. Its default is V/E loading;
passing `PortableLoadMode::VerticesOnly` selects V-only loading, returning an
empty E and communication order. Sampling uses that mode. The graph and E index
remain structurally checked in the manifest; missing or corrupt E payload files
do not block sampling, but do fail default estimator restart. The reported state
digest identifies the source manifest, not a verification of unread E payloads.
`SaveTPSPortable` accepts the
vector TPS and ordered graph; it has no application-metadata argument.
These functions are
collective except descriptor reading and the pure partition helper. Context
and communicator belong to the caller.

The state format supports a fixed simple graph, one physical axis per V, complete
vector TPS, supported real/complex float/double codec and no empty reader ranks.
The checkpoint supports explicit isolated sites in the library; the current
sampler's line/group graph restrictions still apply. CPU gqten native records
require matching integer widths, endianness and scalar type. CUDA records use
the existing TCAPI CUDA format. CPU/CUDA interchange is rejected.

Only the completed format described here is accepted. An unknown common format/convention is rejected.
Additional algorithm state, ancilla marginalization, MCMC-chain
resume, automatic resource selection and dynamic MPI resizing are outside this format.

The CPU acceptance suite exercises independent asymmetric E, noncontinuous
site labels, edge order, 1/2/3/4-rank reads, alternate owners, actual estimator
complex cluster-state output, changed-MPI sampling statistics, restart and
corruption rejection. ROQUO GPU validation passed on GB200 with one shared GPU (1/2/4 MPI ranks)
and two GPUs (1/2 ranks), including loading a checkpoint from a previous job
on another node. All 32 GPU cases passed, including byte-preserving reload,
complex-state sampling and next-cycle resume. Simultaneous multi-node execution
and large-state performance remain unvalidated.


The V-only sampler path was separately validated on ROQUO GB200 with 1/2/4
ranks sharing one GPU (17 cases). Removing the entire message directory left
sampling and its seeded output unchanged; estimator restart still rejected
missing/corrupt E files. Complete estimator restart preserved all V/E bytes.

The CPU suite also checks state-only inputs,
optional estimator resume files, rejection of damaged/mismatched app files by
estimator and their complete independence from sampling. The state/app-file split, changed-circuit estimator restart and custom ownership
were also validated on ROQUO GB200: 83 total cases passed across one shared GPU
and two GPUs, including a later job reading the same state on another node.
Same-circuit steps continued, changed circuits started at zero, bare state input
worked without an app file, corrupt app files affected only estimator, and
V-only sampling remained byte-for-byte identical with E files absent.

The estimator-owned atomic publication wrapper and disconnected-gate rejection
were subsequently checked with 19 fault/topology regression cases on CPU and
on both one shared GPU and two GPUs. The final directory remained absent while
state-only staging was visible; metadata/write/publication failures did not
publish a half-checkpoint, and completed checkpoints resumed at the saved step.
Old app `.pending` input was rejected while genuine library-only input remained
valid. Component-local and single-site gates on disconnected graphs, and routed
gates within a connected graph, passed. This follow-up GPU run passed 73 cases
in total, including the existing changed-MPI restart and V-only sampler checks.

## Path API and compiler dependency

All checkpoint paths (`ReadTPSCheckpoint`, `SaveTPSPortable`,
`LoadTPSPortable`, and the descriptor's `directory`) use `std::string`.
Estimator resume helpers use the same string paths. The TNBP checkpoint layer
itself has no direct `std::filesystem` include or path type.
Directory/stat/rename operations use POSIX helpers in
`tnbp/framework/fileutility.h`; TCAPI tensor I/O still receives filename strings.
The on-disk format and atomic-publication protocol are unchanged.
The implementation still requires the project's existing C++17 features; this
change does not by itself confirm a build with a specific Fugaku compiler.
The pinned CUDA TCAPI backend also avoids `<filesystem>` for string paths and
streams. Direct callers passing `std::filesystem::path` to CUDA TCAPI must include
`<tcapi/input_and_output_filesystem.h>` explicitly; TNBP's string-based calls do
not need that adapter. The CPU min-TCAPI backend retains its existing filesystem
dependency.

## Tensor filenames

Version 3 uses zero-padded site/edge IDs. The minimum width is six digits;
site width is enlarged to fit the largest site ID in that checkpoint, and edge
width is independently enlarged to fit the largest edge ID. All filenames in
each group use the same width, so nonnegative IDs sort numerically with `ls`.
For example, a site ID of 1000000 makes every site filename seven digits wide:
`site-0000001.ten`, ..., `site-1000000.ten`. No six-digit ID limit is imposed.
Logical IDs and graph/axis order are unchanged. No intermediate-format
compatibility reader or migration is provided. MPI rank numbers do not appear
in filenames.

The completed format, string/POSIX paths and adaptive filename widths were
validated on ROQUO GB200 in 92 cases: one shared GPU with 1/2/4 MPI ranks and
two GPUs with 1/2 ranks, including a later job on another node. A two-rank
estimator state was read with an arbitrary four-rank map and preserved every
V/E byte. Sparse labels `{0,1,9,1000000}` and signed labels `{-4,2,7,20}` also
passed changed-MPI V/E roundtrips and custom-owner V-only sampling. The sparse
case verified seven-digit filenames with only four physical sites.
No intermediate-format compatibility or migration is included.
