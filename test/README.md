# TNBP tests

Each numerical suite owns its sources, fixtures, CMake registration and README.

| Directory | Coverage | CTest filter |
| --- | --- | --- |
| [graph_tensor_contraction](graph_tensor_contraction/README.md) | Graph contraction, axis fusion, bond relocation, graph-only plans and invalid inputs | `graph_contraction` |
| [boundary_mps](boundary_mps/README.md) | Ordered boundary sweeps, absorption and independence from the upper region | `bmps` |
| [mcmc_sampling](mcmc_sampling/README.md) | Grouping, environments, Metropolis/heat bath and sampling integration | `mcmc` |
| [mcmc_sampling/heavyhex51](mcmc_sampling/heavyhex51/README.md) | Optional 51-site graph/circuit experiments | explicit executable runs |
| [backend_integration](backend_integration/README.md) | MPI/backend operations, BP, truncation, QASM and header checks | `integration|header_link` |
| [common](common/README.md) | Shared deterministic fixtures and independent reference calculations | no standalone tests |

The pre-existing `truncation`, `square_root_and_inverse` and `min-tcapi`
directories remain in place; they are not newly registered by this reorganization.

## Build and run

From the repository root, configure the common test entry point:

```sh
cmake -S test -B build/tests -DTNBP_TEST_BACKEND=cpu \
  -DGQTEN_ROOT=/path/to/Source/gqten -DHPTT_LIBRARY=/path/to/libhptt.a
cmake --build build/tests -j 4
ctest --test-dir build/tests -j 3 --output-on-failure
```

Initialize dependencies/build HPTT beforehand. Supply platform-specific
OpenMP, BLAS and LAPACK settings where needed; see the
[backend build instructions](backend_integration/README.md#build).
All executables are placed directly in the build directory. Generated suite
files live below `suites/`, avoiding collisions with executable names.

`cmake -S test/backend_integration ...` remains a compatibility entry point and
registers the same suites, so existing build directories can be regenerated.
Use `-S test` for new builds. Backend dependency setup is shared in
`cmake/Backend.cmake`; each suite has its own `CMakeLists.txt`.

CPU CTest uses MPI 1/2 ranks for sampling and 1/2/3 for graph/boundary shifts,
with one OpenMP/BLAS thread per process. Configure `TNBP_TEST_BACKEND=cuda`
to build CUDA versions of the integration, graph, boundary, sampling CLI,
sampling regression and heavyhex51 executables. Launch them explicitly inside
a scheduler allocation with one visible GPU per process (for example,
`srun --mpi=pmix --gpus-per-task=1 --gpu-bind=single:1`). CPU-only CTest
registrations, including checkpoint fixture generation, are not enabled for CUDA.
Graph contraction, boundary shifts, upper-region independence, and real/complex
sampling CLI/checkpoint tests have passed on one GPU and same-node two GPUs.
Small chain/square distributions and manual/automatic heat-bath distributions
were checked with 4096 retained samples. The 51-site manual/automatic GPU cases
use 64 samples, with independent bitstring checks for product states and execution
checks for D=2. This does not cover every CPU CTest mode on CUDA, multi-node runs,
or general sampling accuracy. CUDA checkpoint fixture executables are built;
run each with one process to generate backend-native files before loading them.
