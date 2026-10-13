# Boundary-MPS tests

This directory owns tests of the boundary operations in
`tnbp/ptns/function/bmps.h`.

## Files and coverage

| Source | Executable | Coverage |
| --- | --- | --- |
| `shift_regression.cc` | `bmps_shift` | Right-to-left long-bond relocation on ordered lines |
| `upper_independence.cc` | `bmps_upper_independence` | Lower absorption and BP/truncation with and without an already contracted upper region |
| `main.cc` and Makefile utilities | Makefile driver | Graph-input boundary-MPS experiment |

`shift_regression.cc` checks real/complex tensors on triangles, missing neighbor
bonds, multiple shared neighbors, multihop shifts, scalar remainders, zero
tensors, existing chains, reversed edges, single-site and empty lines. It
compares all open-boundary components with independent dense enumeration,
checks exact ordered graph-only/tensor edges and owner preservation, and
requires a second sweep to leave both tensors and graph unchanged. Invalid
lines, graphs, missing tensors and mismatched dimensions are rejected.
The generic single-bond API and direct contraction tests are covered in
[`graph_tensor_contraction`](../graph_tensor_contraction/README.md).

`upper_independence.cc` uses an eight-site/four-line graph with diagonals.
Only the lower two lines' tensors are supplied for the lower update. It checks
automatic/manual groups and both edge orders before and after BP/truncation to
bond dimension one, with tensor tolerance `1e-11 * (1 + abs(reference))`.

The `main.cc` driver, Makefile, input data and plotting utilities support
standalone graph-input experiments. They are separate from the CTest targets;
their Makefile configuration is outside the validation scope reported here.

## Run

Configure through [the common test entry point](../README.md#build-and-run), then:

```sh
cmake --build build/tests --target bmps_shift bmps_upper_independence
ctest --test-dir build/tests -R bmps --output-on-failure
```

`bmps_shift_1_rank` / `_2_rank` / `_3_rank` use a 90-second timeout.
`bmps_upper_independence_1_rank` / `_2_rank` use 60 seconds.
These CTest registrations use the CPU backend. Both executables also build with
CUDA and have passed explicit runs on one GPU and same-node two GPUs, with one
visible GPU per MPI process. Three-GPU execution has not been validated.

Sampling's end-to-end absorption/environment checks are in
[`mcmc_sampling`](../mcmc_sampling/README.md), because their final assertion is
the sampled distribution rather than an isolated BMPS operation.
