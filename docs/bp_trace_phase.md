# BP trace-phase convention

The default BP algorithm is `trace-phase-gauge-v1`. Each directed rank-2
message is Frobenius-normalized, then multiplied by
`conj(trace(M))/abs(trace(M))` before storage or MPI transfer. Real tensors
use the corresponding sign. The implementation uses `tcapi::trace` and
`tcapi::scale`; `tcapi::get_elem(ctx,traced,{})` reads the logical rank-0
trace value, including any backend lazy scale, using empty coordinates.

The residual is the mean over all directed edges of the Frobenius distance
between unit-norm messages in this same positive-real-trace convention.
It is invariant under a nonzero scalar phase on the messages. The public
vector and `TensorProductState` BP interfaces use the same implementation.
This also covers measurement copies and BP-based boundary-MPS truncation.

This operation preserves norm and singular values. It neither symmetrizes
messages nor projects them onto the positive semidefinite cone. It cannot
repair non-scalar errors or remove approximation errors accumulated during
previous evolution.

A nonpositive/nonfinite norm, nonfinite trace magnitude, or
`abs(trace(M)) <= 64 * epsilon * norm(M)` is an error. Diagnostics identify
the message/edge, norm, trace, threshold and (during BP) MPI rank. Distributed
BP aborts its communicator on a gauge failure to avoid blocked peer receives;
single-rank and standalone helper calls throw. Applications should retain
their step/layer context in their logs.

## Checkpoints and migration

`SaveTPS` writes canonicalized copies of the messages and appends a
`TNBP_BP_METADATA_V1` record after the existing binary tensor payload. It
records algorithm, gauge, residual definition, source algorithm, source
checkpoint identifier and source byte fingerprint. The input state is not
modified. Old TNBP readers stop before the record. CPU and CUDA tensor
payloads remain backend-specific.

`LoadTPS` accepts legacy files without this record, labels them
`legacy-unspecified`, and canonicalizes their messages before returning.
An unknown or malformed metadata record is rejected. The log records both
input and active algorithms and a source fingerprint. These APIs expect
one TPS checkpoint per stream, not concatenated checkpoint records.

For an auditable derived run, pass a `BPCheckpointMetadata` object to
`LoadTPS`, set its `source_checkpoint` to the input path/run identifier,
and pass the same object to `SaveTPS`. `estimator_restart` does this.
For a fresh run set `source_algorithm="fresh"`. Callers that omit metadata
have source provenance marked `unspecified`, rather than an inferred origin.
The stream fingerprint is FNV-1a-64, explicitly non-cryptographic; keep a
SHA-256 manifest of input checkpoint, executable and settings when stronger
integrity guarantees are required. A non-seekable input reports its fingerprint
as `unavailable`.

Changing algorithms on an old checkpoint is a derived run: preserve the
original file and identify its algorithm and hash. Do not interpret the
canonicalization as reproducing the old iteration trajectory or as erasing
past truncation error.

## Regression coverage

`bp_trace_phase` tests real and complex tensors, trace/scale including lazy
scalar values, norm preservation, positive trace, zero/near-zero/nonfinite
trace rejection, scalar-phase invariance of residuals and normalized
observables, a convergent triangle with internal and MPI boundary edges,
legacy/current checkpoint loading, source fingerprints and restart agreement.
The CPU suite runs with one and two MPI ranks; CUDA runs under an explicit
scheduler allocation with one visible GPU per rank. CPU additionally covers
single precision. Use the standard [test build](../test/README.md).
