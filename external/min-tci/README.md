# Minimal gqten TCAPI adapter

The C++17 entry point is `tcapi/tcapi.h`, with the `tcapi` namespace and
`gqten::tensor<T>` for `float`, `double`, `std::complex<float>`, and
`std::complex<double>`. Add this directory's `include` directory and the gqten
include directories to the compiler search path.

This is an incremental migration. The existing `tci/tci.h` implementation remains
available while TNBP callers migrate. Including both headers is supported; their
context handles are distinct. The directory name `min-tci` is retained for now.

The specification baseline is the [TCAPI development specification](https://tensorcomputingapi.github.io/dev/contents/specification/index.html),
retrieved on 2026-10-04, HTML SHA-256
`e432a780c142a69ffa937a37ff1d6fd816d7c98ee9f264f55c3118b1d5132b4d`.
`version<Ten>()` returns "1.0", the user-selected specification baseline, not
a backend release number or a claim of complete support. The limitations below
remain in effect. General `eig` and `eigvals` are deferred.

## Implemented API subset

- Associated types, `List`, `CRef`, `Pair`, and `Map` (`std::unordered_map`).
- `create_context`, `destroy_context`, `version`, `show`.
- `load`, `save` (paths and binary streams).
- `order`, `shape`, `size`, `size_bytes`, `get_elem`.
- `allocate`, `zeros`, `eye`, `fill`, `assign_from_range`, `random`.
- `copy`, `move`, `clear`.
- `set_elem`, `to_range`, `close`, `convert`.
- `reshape`, `transpose`, `cplx_conj`, `to_cplx`, `real`, `imag`.
- Mutable and const `for_each`, `for_each_with_coors`.
- `expand`, `shrink`, `extract_sub`, `replace_sub`, `concatenate`, `stack`.
- `diag`, `norm`, `normalize`, `scale`, `trace`, `linear_combine`, `contract`, `svd`, `trunc_svd`, `qr`, `lq`, `eigh`, `eigvalsh`, `exp`, `inverse`.

All other TCAPI functions, including new-namespace I/O and diagnostic verbosity,
are pending subsequent migration units. The legacy TCI implementation is not a
substitute for the missing public TCAPI functions.

## Backend contracts

Create a `tcapi::context_handle_t<TenT>` and call `tcapi::create_context(ctx)` before
using it. Destroying the context invalidates API use until reinitialization;
CPU tensor storage is owned independently and remains destructible. Creating an
already active CPU context and destroying an inactive one are idempotent.

Shape and coordinate components use gqten's signed 32-bit type. Negative or zero
bond dimensions are rejected with `std::invalid_argument`. Empty shape `{}` is a
valid scalar with logical size 1. Counts and byte products are checked for
`size_t` overflow before allocation; future linear algebra routines must also
check their BLAS/LAPACK integer limits.

`size_bytes` reports numeric element storage, excluding object metadata and
allocator overhead. For scalars this includes one inline element in gqten's
scale field; dense arrays report the heap element payload. gqten's separate
lazy scaling factor is treated as representation metadata for dense arrays.

The default gqten tensor is scalar 1 with no heap payload. `clear` releases the
payload and restores that backend default; `move` transfers ownership and
restores its source to the same default. Portable callers must use `tcapi::copy`
for duplication even though gqten also supports ordinary C++ copies. Returned
tensors can be used to initialize values normally. Direct gqten assignment of a
scalar into a dense tensor can retain the old allocation: call `tcapi::clear`
before such backend assignment if immediate release is needed.

`allocate` provides no initial-value guarantee. Write elements before reading
them. `copy` duplicates the representation, including backend bond labels and
lazy scale, without evaluating uninitialized element values.

Range maps are invoked once per logical element, including once with `{}` for a
scalar. Traversal uses the first bond as the fastest-varying coordinate. Negative
input offsets are rejected, but input range length is not provided by the API;
the caller must supply a valid random-access range. Output map offsets must lie
in `[0, size)`, and the caller must provide storage for that many elements.
`random` stores values produced directly by `gen()`, without introducing a
backend random distribution.

Tests live in `test/min-tcapi/` at the repository root. The gqten submodule is not
modified by this adapter.

## Element manipulation

Dense in-place `reshape` retains the payload and its column-major linear order.
Scalar `{}` can be reshaped to or from any positive shape with product 1.
`transpose` validates a complete permutation before using gqten/HPTT. Output
forms of reshape, transpose and conjugation support an output that aliases the
input, and otherwise produce independent storage. Replacing a dense output with
a scalar releases its old payload. Shape/permutation errors preserve the output.

Complex conjugation includes gqten's lazy scale. Component extraction uses
logical values, including complex scaling; already-complex `to_cplx` and
already-real `real` return deep copies. Real-input `imag` returns zeros.

Both traversal functions return void. Mutable callbacks receive a temporary
logical element by mutable reference; const callbacks receive a const element.
Coordinates are const, and callbacks visit scalars once with empty coordinates.
References supplied to callbacks must not be retained beyond the callback;
callbacks must not independently modify the input tensor through captured aliases.
A mutable traversal stages its output in an additional tensor and commits only
when all callbacks succeed, using O(size) extra storage. If a callback throws,
the tensor is unchanged; external callback side effects cannot be rolled back.
Const traversal does not duplicate or replace the tensor.

## Regions and joining tensors

`expand` appends zero-filled regions; increments must be nonnegative and their
sums must fit `bond_dim_t`. `shrink` uses selected-axis half-open ranges and
`extract_sub` requires one half-open range per axis. Empty slices are rejected
under this backend's no-zero-dimension restriction. `replace_sub` requires equal
input/subtensor order and an origin whose entire replacement region fits.

All region output forms build a temporary before replacing the output. The
output may alias either input, including the subtensor in `replace_sub`.
Validation failures leave output tensors unchanged. Region operations retain
existing backend bond labels. Scalars support empty-map expansion/shrink,
empty-range extraction and scalar replacement with empty coordinates.

`concatenate` and `stack` accept `List<CRef<TenT>>` (for example,
`{std::cref(a), std::cref(b)}` with an explicit tensor template argument).
Inputs remain unchanged and may be referenced more than once. Results own
independent storage and have no backend bond labels. Empty input lists are
rejected. Concatenation requires an existing axis and matching other dimensions;
stack inserts a new axis at any position from 0 through the input order, including
stacking scalars into a vector. Scalar concatenation has no valid axis.

## Comparison and type conversion

`close` requires identical shapes and an absolute difference at most epsilon for
every logical element. Complex differences use their modulus (hypot), with no
relative term or sum over elements. Negative/NaN epsilon is rejected; positive
infinity is allowed. Nonfinite input components compare false, even against
themselves. These nonfinite conventions match the reference CUDA backend.

`convert` supports all 16 ordered pairs of the four gqten element types, using
live CPU source and destination contexts (which may differ). Shape and logical
values are preserved subject to precision conversion. Real-to-complex adds zero
imaginary parts; complex-to-real discards the imaginary part; complex-to-complex
converts both components independently. Equal types use a deep copy, including
when input and output are the same object. Backend labels are retained, and
validation/allocation failures before replacement preserve the old output.

Component casts use CPU floating-point conversion. Values beyond the finite
range of the destination type are explicitly mapped to signed infinity. NaN,
infinity and signed zero are supported, without promises about NaN payloads or
identical CPU/GPU rounding. Unit lazy scales are bypassed when fetching numeric
values, avoiding complex multiplication contamination of nonfinite components;
non-unit scales follow gqten arithmetic. CPU/GPU transfers are not implemented.

## Diagonals and scalar arithmetic

`diag` converts a vector to a dense diagonal matrix, or extracts the main
diagonal of a rectangular matrix (length min(rows, columns)). Other ranks are
rejected. Logical scaling is included; off-diagonal elements are zero. Both
forms allocate an independent result before replacement and permit output aliasing.

`norm` reuses gqten's BLAS-backed Frobenius norm, including its lazy scale and
scalar absolute value. Dense vector length must fit the BLAS `int` interface.
Extreme-value and nonfinite results follow gqten/BLAS floating-point arithmetic,
including rounding/overflow in the product of raw norm and absolute lazy scale.

`scale` reuses gqten's lazy factor multiplication. Separate outputs are deep
copies; an identical input/output uses the in-place operation. `normalize`
returns the original norm and usually divides the lazy factor by it. If that
factor division overflows or underflows to zero, it instead normalizes logical
elements into a temporary result. Zero, infinite and NaN norms raise
`std::domain_error` before mutation. Invalid output-form calls preserve the
existing output. Normalization is subject to floating-point rounding, especially
for subnormal values. The scalar/default-state and zero-dimension policies above
continue to apply.

## Trace and linear combinations

`trace` accepts disjoint pairs of distinct, equal-dimension axes. Pair direction
and pair ordering are unrestricted; remaining axes retain their original order.
Full trace returns a scalar, and an empty pair list returns an independent copy.
Output forms permit aliasing and preserve the output on validation failure.
The adapter follows gqten's coordinate-sum algorithm rather than including its
trace header, which contains a non-inline free helper incompatible with this
public header's multiple-translation-unit use. No gqten source is modified.

`linear_combine` accepts non-owning const references:

```cpp
auto sum = tcapi::linear_combine<Ten>(ctx, {std::cref(a), std::cref(b)});
auto difference = tcapi::linear_combine<Ten>(
    ctx, {std::cref(a), std::cref(b)}, {Elem(1), Elem(-1)});
```

Coefficient omission means all ones. Inputs must be nonempty and identically
shaped; explicit coefficient count must match. Repeated references are allowed.
Inputs remain unchanged and the returned result owns independent storage.
Dense combinations reuse gqten::LinearCombine and its BLAS int-length limit;
scalar combinations are handled by the adapter because gqten's dense routine
assumes allocated array data. Internal pointer adaptation matches gqten's old
non-const input signature, whose implementation only reads the inputs. Numerical
rounding and nonfinite arithmetic follow the underlying CPU operations.

## Contraction

`contract` keeps the legacy label remapping and gqten::Contract kernel. Inputs
are const tensor references (not CRef lists). Integer labels may be any int32
value; string_view labels mean one label per byte/character, including spaces
and punctuation, without the legacy comma/space token parsing extension.
Repeated labels within an operand are rejected. Every free label must appear
once in the output; its position defines the output axis order. Shared input
labels must have equal dimensions and be contracted: retaining a shared label
in the output is unsupported by this backend and raises invalid_argument.
Contraction is bilinear, without implicit complex conjugation.

A temporary result permits output aliasing either or both inputs. Scalar inputs
use an adapter multiplication/permutation because the native routine requires
dense input storage. Full contraction produces a scalar. Validation failures
leave output unchanged. Input and output element counts are limited to INT_MAX
because native matrix dimension/allocation products use signed int. Floating
point rounding and nonfinite arithmetic follow gqten/BLAS. The gqten submodule
is unchanged.

## Singular value decomposition

`svd` and the two specification overloads of `trunc_svd` use the first `rows`
axes as matrix rows, requiring `1 <= rows < rank`. U and V-dagger retain the
corresponding input axes, with a shared decomposition axis appended/prepended.
Sigma is a **dense real rank-2 diagonal matrix**, not a vector. Its storage is
O(chi squared). Outputs must be distinct from each other; an output may alias
the input. All factors are prepared before outputs are replaced.

The adapter uses gqten::TruncSVD with all singular values retained to obtain
LAPACK error reporting (gqten::SVD discards its status). The explicit `svd`
(gesvd) driver avoids the native `auto` fallback's incomplete byte copy.
The gqten submodule remains unchanged. Input element counts must fit signed int;
nonfinite inputs are rejected, and reported numerical failures become exceptions.

Truncation selection is implemented in the adapter to follow the specification:
values below s_min are discarded, chi_min never restores them, and the smallest
allowed rank meeting the inclusive target-error bound is kept, up to chi_max.
The short overload is equivalent to chi_min=1 and target=0. Parameters require
1 <= chi_min <= chi_max, finite nonnegative s_min and target. No surviving
singular value raises domain_error under the existing zero-dimension restriction.
Squared weights are normalized by the largest singular value and accumulated
as long-double suffix sums to avoid overflow and cancellation. For the zero
matrix, relative truncation error is defined as zero; s_min=0 keeps at least
chi_min values where available. Validation/numerical failures preserve outputs
and trunc_err. Full factors are computed before truncation, so peak storage is
larger than the retained factors alone. The legacy s_min-only overload is not
part of the new public API.

## QR and LQ

`qr` reuses gqten::QR. `lq` follows the legacy implementation: transpose the
row/column axis groups, apply QR, then transpose the two factors back. Ordinary
transpose suffices for complex inputs as well: A^T=Q'R' gives A=R'^T Q'^T,
with orthonormal rows in Q'^T. No conjugation or phase convention is imposed.
Both return thin factors with shared dimension min(product(row axes),
product(column axes)); axes within each input group retain their order.

Both require 1 <= row bonds < rank, finite input values, and an input element
count fitting signed int. Outputs must be distinct, but either may alias the
input. Temporary factors preserve existing outputs on validation failure and
avoid aliasing native buffers. LQ also copies/transposes the input. Numerical
LAPACK status handling remains gqten's assert-based behavior; the native QR API
does not expose status codes. No gqten, legacy tci, or TNBP code is changed.

## Symmetric / Hermitian eigenproblems

`eigh` and `eigvalsh` require a real symmetric or complex Hermitian input after
matricization. General `eig` and `eigvals` are deferred. The row-axis count must
satisfy 1 <= rows < rank and the two axis groups must have equal products.
Symmetry/Hermiticity is a caller precondition, not numerically checked or repaired;
the native solver uses the upper triangle (LAPACK UPLO=U).

`eigh` returns ascending eigenvalues in a dense real rank-2 diagonal matrix
lambda_mat of shape {n,n}, and eigenvectors of shape {row axes...,n}.
`eigvalsh` returns an ascending real rank-1 list of shape {n}, requesting no
eigenvectors from LAPACK. The native routine still uses an n-by-n work buffer.
Inputs are materialized with unit scale before gqten::EigHerm so negative lazy
scales preserve ascending order and complex storage factors are interpreted via
the logical Hermitian matrix. Finite inputs and signed-int element-count limits
are checked. Output aliasing the input is supported; eigh outputs must be distinct.
Validation failures leave outputs unchanged. Repeated eigenvalues have no fixed
basis/phase convention. Native LAPACK status handling remains assert-based, as
with QR; EigHerm does not return a numerical status. gqten is unchanged.

## Matrix exponential

`exp` supports real symmetric / complex Hermitian matrices only, with the same
matricization validation and caller precondition as `eigh`. Both in-place and
output forms preserve the original tensor shape. The existing gqten::ExpHermExact
algorithm diagonalizes the logical input, forms B=V exp(lambda/2), and returns
B B-dagger. A temporary output supports aliasing and preserves the output on
validation failure. This is a matrix exponential, not an elementwise operation.
Native eigensolver status handling remains assert-based. Floating-point
exponential overflow/underflow and BLAS propagation are inherited from gqten;
no general non-Hermitian algorithm or overflow guarantee is provided.

## Matrix inverse

`inverse` supports general invertible real/complex square matricizations, with
1 <= row bonds < rank. It refolds the inverse into the original tensor shape.
Both forms support input/output aliasing. Finite logical input values (including
lazy scale) are copied into a work buffer; outputs are replaced only on success.
The adapter calls the same LAPACK GETRF/GETRI routines used by gqten, checking
their status directly because gqten::Inverse does not return it. A singular
matrix raises domain_error; rejected LAPACK arguments raise runtime_error.
Workspace query results are checked before conversion/allocation. Input element
counts must fit signed int. No gqten source is changed.

No tolerance-based singularity cutoff or conditioning estimate is applied:
a nonzero pivot is accepted, and ill-conditioned matrices can lose accuracy.
Overflow/underflow in numerical results is not separately trapped. Validation
or reported LAPACK failure leaves input/output unchanged, including in-place use.

## Binary I/O

`load<Ten>(ctx, storage)` and `save(ctx, a, storage)` accept std::string,
std::string_view (without requiring a trailing null), const char*, and
std::filesystem::path. They also support std::istream/std::ostream and derived
streams such as ifstream/ofstream/stringstream. Open file streams in binary mode.
Each stream call reads/writes exactly one tensor from the current position;
it does not close, seek, or explicitly flush the caller's stream, or consume
following records. This supports consecutive tensor dumps. A failed operation
can advance the stream and leave a partial record; no stream rollback is promised.
Path save truncates/replaces the file, is not atomic, and checks close failures.

The format remains gqten's native sequence: size_t rank, int32 shape entries,
raw element-type scale, then raw elements (none for a scalar). Bond labels are
not serialized. It has no magic/version/type tag and uses host byte order and
native size_t/element representations. Loading requires the same element type
and compatible ABI as the writer; type mismatch cannot reliably be detected.
This preserves legacy gqten/CPU files, not cross-backend or CUDA format compatibility.

The new reader validates dimensions, counts and read completion instead of
using unchecked gqten::StreamRead. Truncation, stream failure and path-open errors
raise ios_base::failure; invalid dimensions and size overflow use the adapter's
usual exceptions. Resource exhaustion can still raise bad_alloc. The writer
retains gqten::StreamWrite and checks stream state. Caller exception masks are
preserved. No gqten or CUDA implementation is modified.

## Display and diagnostics

`show(ctx, a)` prints gqten's human-readable tensor display to standard output,
including logical scaled values (with gqten's default display truncation).
`version<Ten>()` returns "1.0" for the supported tensor types.

TCAPI_VERBOSE follows the reference CUDA implementation: unset, 0 or invalid
values disable adapter logging; 1 emits a single line per public call to stderr;
2 adds elapsed_ms measured with steady_clock. The environment is read once on
first use per process. Lines contain the function name, selected input metadata
and status=ok/exception. Timing excludes metadata formatting and log output;
it includes the operation and cleanup, and is not a benchmark guarantee.

All 53 implemented function names are instrumented. Internal nested calls are
suppressed; calls from user callbacks are separate entries. Tensor metadata is
captured before mutation without reading elements. Logging is best-effort and
cannot replace computational exceptions. A shared mutex prevents adapter log
lines from interleaving across threads. This does not make tensor mutations
thread-safe. Existing gqten messages/assertions are outside TCAPI_VERBOSE control.
The logging implementation is adapted from tcapi-cuda's detail/verbose.h.
