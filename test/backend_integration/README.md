# CPU/CUDA TNBP integration checks

Small deterministic checks of MPO application to MPS, belief propagation and
truncation. Uses `complex<double>` through public TCAPI calls; the same source
selects gqten or tcapi-cuda at build time. Sampling and boundary MPS are excluded.

Coverage:

- Existing `MpiBcast`, `MpiSend` and `MpiRecv`: asymmetric 2x3 tensor, vector and
  scalar; complex values, shape preservation and conjugation after reception.
- One QASM2 Bell circuit: parsing, edge extraction and both `QasmToTPO` overloads.
  `AbsorbTPO` is checked against all four exact state amplitudes. The grouped
  overload also exercises `OptTPObySVD` and rank-local `AttachTPO`.
- `BeliefPropagation` and `BeliefPropagationCondition` on a two-site tree;
  `Truncation` preserves the normalized Bell state with maximum bond dimension 2.
- Nondegenerate Schmidt weights 2:1: BP followed by rank-one truncation gives
  discarded weight 0.2 and normalized state |00>. This also exercises real
  spectrum transfer inside the complex calculation on two ranks.
- One unit-coefficient Pauli string ZX: qubit ordering, shape and all 16 dense
  entries. Coefficient handling, identity-only terms, Y, and file parsing are
  not covered by this smoke test.

State comparisons reconstruct the two-site state by an independent host sum.
After truncation they remove normalization/global phase, avoiding comparison of
arbitrary SVD gauges. Pre-truncation Bell amplitudes retain an absolute check.
Failures abort the MPI communicator. CPU CTest runs have a 120-second timeout.
These are small integration fixtures, not exhaustive numerical or performance
coverage; multi-edge BP convergence and all parser gates remain untested.

## Build

Backend repositories require access permission; see [repository access](../../README.md#backend-repository-access).

CPU (initialize dependencies and build HPTT separately):

```sh
cmake -S test/backend_integration -B build/backend-cpu \
  -DTNBP_TEST_BACKEND=cpu -DGQTEN_ROOT=/path/to/Source/gqten \
  -DHPTT_LIBRARY=/path/to/libhptt.a
cmake --build build/backend-cpu
ctest --test-dir build/backend-cpu --output-on-failure
```

Supply platform-specific OpenMP/BLAS/LAPACK settings to CMake if needed.

CUDA (build on the allocated GPU compute node):

```sh
cmake -S test/backend_integration -B build/backend-cuda \
  -DTNBP_TEST_BACKEND=cuda \
  -DCUTENSOR_INCLUDE_DIR=/path/to/cutensor/include \
  -DCUTENSOR_LIBRARY=/path/to/libcutensor.so
cmake --build build/backend-cuda
```

The CUDA build excludes the CPU adapter's include directory. It uses C++17 and
links MPI, CUDA runtime, cuBLAS, cuSOLVER and cuTENSOR. Initialize the pinned
`external/tcapi-cuda` submodule before configuring. `TCAPI_CUDA_ROOT` defaults
to that submodule; override it to use another checkout and record its commit.

CUDA launch must expose exactly one GPU per rank. For two ranks, the test checks
that both run on the same host and have different GPU UUIDs. Both processes can
legitimately call their own device "0". It records UUID and CUDA_VISIBLE_DEVICES
and fails instead of silently sharing one device. A representative Slurm step
inside an appropriate allocation is:

```sh
srun --nodes=1 --ntasks=2 --gpus-per-task=1 --gpu-bind=single:1 \
  ./build/backend-cuda/backend_integration
```

Scheduler allocation, MPI launch/module configuration and time limit must be
adapted to the site. The generic command above omits site-specific settings.
On 2026-10-05, ROQUO runs with explicit `--mpi=pmix --cpus-per-task=36`
passed: one rank/one GB200 (67 checks) and two ranks/two GB200s (72 checks per
rank), with distinct UUIDs on the same host. Toolchain: GCC14.2.1, CUDA13.3.73,
cuTENSOR2.6.0 and HPC-X2.50/OpenMPI5, tcapi-cuda commit `0cc6394`.
CPU runs passed 66 checks on one rank and 69 per rank on two ranks; extra CUDA
checks validate device isolation. No CUDA-aware MPI is required: TNBP's existing
transfer helpers stage tensor data through host vectors. This does not establish
multi-node correctness, performance, or exhaustive API coverage.

## Header independence

The default build also builds `header_checks`: one translation unit per public
`include/tnbp/**/*.h`, `include/qasm/**/*.h`, and `include/pauli/*.h` header
(37 total), with that header as its first and only include,
plus a translation unit including all headers in reverse lexical order. This
covers individual headers and `inc_all.h` entry points with either backend.
No forced includes or precompiled headers are used. The target compiles object
files without linking them; it checks include dependencies, not cross-unit ODR
correctness or every template instantiation. It does not run sampling routines.

Each header includes its own required standard headers and internal
prerequisites. `parser/qasmutility.h` contains shared `OpQubitCount` metadata so
`qasmtoedge.h` does not require the tensor converter or a tensor backend first.

## Header-only linkage

`header_link_check` compiles three separate sources that include `tnbp/tnbp.h`
and links them into one executable. On GCC/Clang it uses `-O0 -fno-inline` and
explicitly disables unity builds. The test compares the addresses of11 public
non-template functions across source files, ensuring a single shared entity
rather than hiding duplicate definitions with `static`. It also checks a small
get_range example. This adds no production library or library-build requirement.

Before the inline fix, the test failed to link with11 duplicate symbols:
get_range, string MpiBcast, seven lattice-generation functions, SitesFromQasm and
EdgesFromQasm. These header-defined functions now have external inline linkage.
Template functions and in-class definitions already have appropriate header-only
semantics; explicit MPI datatype specializations were already inline.
