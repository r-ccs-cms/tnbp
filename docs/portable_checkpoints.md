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

Estimator step information is optional for TPS input. Sampling does not require
estimator circuit/settings metadata or adopt its truncation parameters. It
prints the source state digest and, when available, estimator next-step label.
Portable sampling seeds independent rank streams from `(seed, rank)`; changing
MPI count therefore changes the random sequence, while retaining the target
state. Legacy sampling retains its previous seed behavior.

To continue estimator evolution instead:

```sh
mpirun -np 4 ./estimator --circuit circuit.qasm --sparse_pauli observables.txt \
  --portable_load /shared/run/checkpoint-10 --step_end 20 \
  --portable_save /shared/run/checkpoint-20
```

Estimator validates its application version, complete-cycle phase, circuit
content SHA256, BP/truncation settings and next step. An explicit `--step_start`
must equal the saved next step; otherwise it is inferred. Set `--step_end` to the
exclusive endpoint. Loading with `step_end == next_step` saves without evolution.
Legacy rank-local `--loadname`/`--savename` continue to use their existing format;
legacy-to-portable conversion is not implemented here.

## Layout and integrity

```text
checkpoint-10/
  manifest.txt
  COMMITTED
  vertices/site-0.ten
  vertices/site-1.ten
  ...
  messages/edge-0-forward.ten
  messages/edge-0-reverse.ten
  ...
```

`manifest.txt` starts with `TNBP_PORTABLE_TPS 1`. It contains the codec/ABI,
BP conventions, axis convention, sorted site IDs, ordered edges with stable
snapshot edge IDs, file records, writer-count provenance and optional quoted
application key/value fields. File records contain path, length, SHA256 and
shape. All counts are checked decimal integers; strings use C++ `std::quoted`.
No JSON, external parser, hash command, or new library dependency is used.
SHA256 is implemented incrementally in the library and checked against known
vectors by the regression suite.

Vertex virtual axes follow the saved incident-edge order; the single physical
axis is last. For saved edge `(a,b)`, `forward` means `a` to `b`, `reverse` means
`b` to `a`. E axis zero is ket, axis one is bra. Directions are independent
records. The loader never sorts/reverses graph edges or transforms tensors.
The state digest covers codec, conventions, graph and logical record hashes;
it excludes MPI ownership, writer count and application metadata.

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

`ReadTPSCheckpoint(directory)` returns the graph, logical records and optional
application fields. `PortableContiguousOwners(sites,comm)` supplies the CLI's
balanced contiguous partition. `LoadTPSPortable<T>(ctx,directory,owners,comm)`
also accepts an explicit owner map and returns local V/E, SiteIdx/EdgeIdx,
new owners, communication order and descriptor. Its default is V/E loading;
passing `PortableLoadMode::VerticesOnly` selects V-only loading, returning an
empty E and communication order. Sampling uses that mode. The graph and E index
remain structurally checked in the manifest; missing or corrupt E payload files
do not block sampling, but do fail default estimator restart. The reported state
digest identifies the source manifest, not a verification of unread E payloads.
`SaveTPSPortable` accepts the
vector TPS, ordered graph and optional application fields. These functions are
collective except descriptor reading and the pure partition helper. Context
and communicator belong to the caller.

Version 1 supports a fixed simple graph, one physical axis per V, complete
vector TPS, supported real/complex float/double codec and no empty reader ranks.
The checkpoint supports explicit isolated sites in the library; the current
sampler's line/group graph restrictions still apply. CPU gqten native records
require matching integer widths, endianness and scalar type. CUDA records use
the existing TCAPI CUDA format. CPU/CUDA interchange is rejected. Unknown
application fields are optional, but an unknown common format/convention is
rejected. Additional algorithm state, ancilla marginalization, MCMC-chain
resume, automatic resource selection and dynamic MPI resizing are outside v1.

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
