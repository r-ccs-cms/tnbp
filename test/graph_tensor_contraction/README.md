# Graph tensor contraction tests

This directory owns tests of `tnbp/framework/graph_contraction.h`.

## Files and coverage

- `regression.cc`: deterministic regression executable. Tests ordered
  contraction plans against explicit expected edges; real/complex contraction
  against independent dense open-boundary values; multiple common neighbors,
  unequal dimensions, reversed edge order, fresh/in-place result labels and
  partial external tensor/owner inputs. Fifteen invalid-contraction cases check
  collective rejection without changes to edges, IDs, owners, shapes or values.
- Generic single-bond relocation tests cover absent destination bonds, shared
  endpoints and scalar remainders. They check every open-boundary component,
  exact tensor/graph-only edge order, owner preservation and unselected edges.
- `main.cc`, `Makefile`, `Configuration`, `run.sh`, input and plotting utilities:
  standalone graph-input contraction driver. It is separate from the CTest
  regression target; its Makefile configuration is outside the validation
  scope reported here.
- `CMakeLists.txt`: this suite's target and MPI test registrations.

The independent reference fixture is shared through
[`../common/graph_fixture.h`](../common/graph_fixture.h). It does not use the
contraction or SVD being tested to compute expected boundary components.
Boundary-line selection and multihop sweep tests belong to
[`boundary_mps`](../boundary_mps/README.md).

## Run

Configure through [the common test entry point](../README.md#build-and-run), then:

```sh
cmake --build build/tests --target graph_contraction_regression
ctest --test-dir build/tests -R graph_contraction --output-on-failure
```

CTest runs `graph_contraction_regression_1_rank`, `_2_rank` and `_3_rank`.
Direct example: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpirun -np 2
build/tests/graph_contraction_regression`.

Open-boundary comparisons use tolerance `1e-10 * (1 + reference scale)`.
Both real and complex tensors are exercised. The executable supports the
CPU/CUDA backend selection. CUDA execution has also passed on one GPU and
same-node two GPUs, with one visible GPU per MPI process. The three-rank CTest
coverage is CPU-only; CUDA executables are launched explicitly in a GPU allocation.
