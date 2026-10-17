# CUDA backends for estimator applications

`estimator` and `estimator_restart` support `USE_TCAPI_CUDA` to select
`tcapi::cuda::Tensor<std::complex<double>>`. `USE_SINGLE` selects
`std::complex<float>` instead. Defining both `USE_TCAPI_CUDA` and `USE_CYTNX`
is rejected. Without either selector the existing gqten backend remains the
default. The CLI `--backend` option still selects graph/device topology; it is
not the tensor backend selector.

The tcapi-cuda and gqten repositories are not publicly available as of October
2026. You need access permission for the backend you choose; the submodule
entry alone does not make its contents available. See [repository access](../README.md#backend-repository-access).
With tcapi-cuda access, initialize its pinned submodule from the TNBP root:

```sh
git submodule update --init --recursive external/tcapi-cuda
```

Build from either application directory, for example:

```sh
cd apps/estimator
make CONFIGURATION=../Configuration.cuda \
  CUDA_ROOT=/usr/local/cuda \
  CUTENSOR_ROOT=/path/to/cutensor
```

`TCAPI_CUDA_ROOT` defaults to `../../external/tcapi-cuda` from either app
directory. Override it to use another checkout.

The same command works from `apps/estimator_restart`. Both existing Makefiles
produce an executable named `estimator` in their own directory. Set
`CUDA_LIBDIR` or `CUTENSOR_LIBDIR` if libraries are not in `lib64` or `lib`,
respectively. `CCCOM` defaults to `mpicxx`. Add `USE_SINGLE=1` for the
single-precision selection. The initial runtime validation uses double precision.

The CUDA configuration includes only tcapi-cuda, CUDA and cuTENSOR alongside
TNBP; it does not include min-tci/gqten. It links cudart, cuBLAS, cuSOLVER and
cuTENSOR, with library search paths also recorded as runtime paths. Use the MPI
compiler associated with the selected CUDA-compatible execution environment.
The ordinary `make` path retains the existing CPU Configuration.

## Standard default-pool backend

The pinned `tcapi-cuda` uses the CUDA device default pool for Tensor storage
and internal workspaces. CUDA memory-pool support is required; there is no
fallback to direct allocation. Only this CUDA implementation is bundled.

`TCAPI_CUDA_VARIANT` has been removed. Old `direct` and `pool` selections are
rejected with migration guidance. Build the standard backend without that
variable. An explicit `TCAPI_CUDA_ROOT` can select another checkout, but such
implementations are outside regular TNBP validation. All implementations share
header/type names: include exactly one and rebuild after switching. The Make
entry points always rebuild their application object before linking.

For CMake tests and the MCMC CLI, use a fresh build directory when changing
backends or checkout paths. `TCAPI_CUDA_ROOT` defaults to
`external/tcapi-cuda`; an existing CMake cache keeps its explicit root until
reconfigured. For example, from the repository root:

```sh
cmake -S test -B build-cuda -DTNBP_TEST_BACKEND=cuda \
  -DCUDAToolkit_ROOT=/usr/local/cuda \
  -DCUTENSOR_INCLUDE_DIR=/path/to/cutensor/include \
  -DCUTENSOR_LIBRARY=/path/to/cutensor/lib/libcutensor.so
cmake --build build-cuda --target mcmc_sampling
```

`create_context(ctx)` initializes execution resources and is required before
operations. The CUDA Context releases them automatically at scope exit;
ordinary application tests need no explicit `destroy_context(ctx)`. Tensor
storage is independently owned and may outlive its creating Context. The entry
points declare Context before tensor containers, so tensors are released first.

Normal data-producing APIs complete numerical work before returning, including
host staging used for blocking MPI transfers. Internal workspace returns may
still be pending. Tensor destruction is not a synchronization barrier for
outstanding GPU uses. External CUDA work needs explicit ordering and completion
before storage is released; metadata queries are not general barriers.

`TCAPI_CUDA_DEFAULT_POOL_RELEASE_THRESHOLD_BYTES` sets the shared default-pool
release threshold before first Context creation. Unset leaves the existing
CUDA setting unchanged; it does not select the maximum threshold. `0` requests
release of unused memory at synchronization, while `18446744073709551615`
requests maximum retention. This is neither preallocation nor an allocation cap,
and affects other users of the same default pool. The dedicated implementation's
`TCAPI_CUDA_POOL_RELEASE_THRESHOLD_BYTES` is ignored by the standard backend.
See [pool configuration](../external/tcapi-cuda/docs/pool-foundation.md).

## GPU assignment

Launch with **one visible GPU per MPI process**. Each application checks this
condition and selects CUDA device0 before creating any TCAPI context. If no GPU
or several GPUs are visible, it stops with an explanatory error. The launcher
is responsible for assigning different physical GPUs to different local ranks.

A representative Slurm step inside an allocation for two GPUs is:

```sh
srun --mpi=pmix --ntasks=2 --gpus-per-task=1 --gpu-bind=single:1 \
  ./estimator --circuit circuit.qasm --sparse_pauli observables.txt
```

Adapt account, CPU resources and MPI settings to the site. On ROQUO the tested
setup uses gcc/14, cuda/13.3, hpcx/2.50 and an explicit
`--cpus-per-task=36` on the step. Both ranks may see a device numbered0 while
using different physical GPUs. MPI tensor traffic still uses host staging;
CUDA-aware MPI is not required.

## Restart scope

`estimator_restart` uses the backend's stream I/O. tcapi-cuda commit `0cc6394`
includes this functionality. Checkpoints should be loaded with the same backend,
precision and rank layout; cross-backend checkpoint conversion is not provided
by this change. Computational options and existing checkpoint naming are
unchanged. Sampling/boundary-MPS functionality is outside this app integration.

## Validation

On October 17–18, 2026, the standard default-pool pin (`9dea1a0`) was checked
in double precision on ROQUO with one GB200 and two GB200s on one node.
Backend/MPI integration, graph contraction, BMPS shifts and upper-region
independence, BP trace phase, and real/complex MCMC CLI/checkpoints passed.
The seven MCMC regression modes also passed, using heat bath for fixed-site
statistical checks. Long cases were split without reducing sample counts or
relaxing tolerances.

All 488 bundled-input expectation values agreed with the CPU reference within
6.2e-14 on both GPU layouts; all four BP layers converged. Estimator and restart
Bell checks passed. Restart preserved the tensor/graph payload on load/resave;
provenance metadata was updated as intended. Saving after one circuit step and
resuming for another matched two uninterrupted steps for all 15 nonidentity
two-qubit Pauli observables, and agreed with an independent dense calculation
within 1e-15. Integration also checks empty measurement requests and multiple
observables on the same graph edge.

These are correctness checks, not a performance comparison, multi-node test or
GPU memory-sanitizer result.

### Historical backend results

The results below used the stated historical backend commits.

On 2026-10-05 both applications built and ran in double precision on ROQUO
with one rank/one GB200 and two ranks/two GB200s on one node. A two-qubit Bell
circuit produced Z=0 and ZZ=1 in both applications. `estimator_restart` also
preserved each rank's checkpoint bytes through load/resave without further
evolution. CPU/gqten two-rank regression passed the same checks. These are small
smoke tests, not coverage of every circuit/option or cross-backend checkpoints.

The bundled 156-qubit kicked Ising input also completed in `estimator` on CPU
with 1/2/4 MPI ranks and on ROQUO with 1/2 GPUs. All 488 expectation values
agreed with the CPU one-rank result to within 8e-14, and all four BP layers
converged. The maximum retained bond dimension was only 5; this is a correctness
smoke test, not a GPU performance benchmark.

On 2026-10-17, the same 156-qubit input was run four times per backend on one
GB200, alternating direct/pool order. The original backend (`0cc6394`) and pool
backend (`507b210`) produced identical printed values for all 488 observables
in all eight runs; all four BP layers converged with identical iteration counts.
Median elapsed time was 39.71 s (direct) versus 24.75 s (pool), a 1.60x speedup.
Timing includes process launch, initialization and teardown, but excludes builds.
This small-bond result does not establish large-bond or multi-GPU performance;
pool runs of the restart and MCMC applications were not included.
