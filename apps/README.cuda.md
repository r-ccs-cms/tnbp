# tcapi-cuda for estimator applications

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
