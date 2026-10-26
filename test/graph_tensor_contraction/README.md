# Graph tensor contraction tests

This directory owns tests of `tnbp/framework/graph_contraction.h`.

## Files and coverage

- `regression.cc`: deterministic regression executable. Tests ordered
  contraction plans against explicit expected edges; real/complex contraction
  against independent dense open-boundary values; multiple common neighbors,
  unequal dimensions, reversed edge order, fresh/in-place result labels and
  partial external tensor/owner inputs. A three-rank case verifies operand owners
  finish while an uninvolved rank waits outside the contraction. Fifteen invalid
  cases run as separate bounded MPI jobs and must abort with a diagnostic.
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

## Shape-only cost evaluation

`graph_contraction_cost_shapes` tests dense MACs and output sizes against
analytic chain/triangle/K4/square cases, all 720 K4 edge-schedule permutations,
common-neighbor shapes, priority/label rules, partial/empty schedules, and
validation/overflow without input mutation. This target has no backend or MPI
link dependency and can also be compiled directly with a C++17 compiler and
`-I ../../include` from this directory.

`graph_contraction_cost_tensor` compares every predicted graph and site shape
with actual binary tensor contractions, then independently runs the numerical
network executor. It covers triangle/K4/square/honeycomb, reversed edges,
priorities/fresh labels, split owners and an empty rank zero. CTest runs it on
CPU MPI 1/2/3 ranks. CUDA execution is not part of this addition's validation.
See [the cost API](../../docs/graph_contraction_cost.md) for inputs and reports.

## Small region plans

`graph_region_plan` compares the planner with an independent leaf-subset oracle
for 1–5 leaves and several repetition counts. It checks weighted common work,
owner normalization, deterministic ties, six-leaf fallback, overflow recovery,
and unchanged outputs on failure. See [the region API](../../docs/graph_region_contraction_plan.md).

`graph_contraction_cost_symbolic` compares 384 small D/chi instances and ordered
K4 contraction schedules with the numerical cost evaluator, including intermediate
common-neighbor shapes and edge order. It also checks invalid input, immutable
inputs, target-owner retention, external-neighbor footprints and model-plan
selection against all prepared candidates. No MPI/backend is needed.

Symbolic region tests compare all 18 four-leaf candidate plans with standalone
schedule evaluation for multiplicities 1, 2 and 5, including costs, result sizes,
owner footprints and resulting edges. Invalid extra shapes must fail without
modifying the caller input.
