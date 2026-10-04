# min-tcapi CPU tests

Tests follow the TCAPI type-system and function categories. Each executable tests
float, double, complex float, and complex double, and exits nonzero on failure.
Runtime checks remain enabled with `NDEBUG`.

| Directory | Current coverage |
|---|---|
| type_system | Associated/auxiliary types, context lifetime, legacy-header coexistence, multiple translation units |
| read_only_queries | Scalar and dense shape/size/bytes, logical values with lazy scaling, invalid coordinates |
| construction_and_destruction | Constructors, direct generator values/counts, independent copies and backend metadata, ownership transfer, actual clear release, invalid shapes/overflow, throwing generator |
| tensor_manipulation | set_elem, scalar/dense reshape, 3D transpose/inverse, output aliasing, conjugation, component extraction, mutable/const traversals, callback exceptions, region padding/slicing/replacement, alias safety, CRef joins and scalar stack |
| miscellaneous | Mapped range construction/export, scalar callbacks, complex round trip, invalid offsets/context; absolute close boundaries/nonfinite inputs; all 16 conversion pairs, deep copy, alias, signed zero and narrowing overflow |

`input_and_output` will be added with its implementation unit. The foundation, tensor-manipulation and initial linear-algebra APIs are currently tested.

The `tensor_linear_algebra` executable covers bidirectional diag (including rectangular matrices), norm/normalize/scale, scalar and complex values, output aliasing and independence, invalid ranks/contexts, zero/nonfinite norm rejection and tiny positive norms. It also covers partial/full/empty trace, reversed/disjoint pairs, aliasing, and CRef linear combinations with omitted/complex coefficients, repeated references, scalar inputs and validation failures. Numerical comparisons use 64 times the real-type epsilon with a mixed absolute/relative bound.

## Build and run

Use `make check` with a C++17 compiler and the required gqten include/link flags.
No MPI initialization or launcher is needed for these foundation tests. The MPI
compiler wrapper supplies headers pulled in by gqten. The Makefile accepts
`CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, `LDLIBS`, and `BUILD_DIR`.

Transpose tests require a built gqten HPTT library. For Apple Clang and Homebrew,
build it from the repository root:

```sh
make -C external/tensor-ng-dev/Source/gqten/ext/hptt -j4 \
  CXX=/usr/bin/clang++ \
  CXX_FLAGS='-std=c++17 -O3 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include'
```

Then, from `test/min-tcapi/`:

```sh
OMPI_CXX=/usr/bin/clang++ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 make -j4 check \
  CXX=/opt/homebrew/opt/open-mpi/bin/mpicxx \
  CXXFLAGS='-std=c++17 -O1 -g -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/openblas/include' \
  LDLIBS='-L../../external/tensor-ng-dev/Source/gqten/ext/hptt/lib -lhptt -L/opt/homebrew/opt/libomp/lib -lomp -L/opt/homebrew/opt/openblas/lib -lopenblas'
```

For memory/undefined-behavior checks use a separate `BUILD_DIR=build-asan` and
append `-fsanitize=address,undefined -fno-omit-frame-pointer` to `CXXFLAGS`.
For release checks use `BUILD_DIR=build-release` and replace `-O1 -g` with
`-O3 -DNDEBUG`. Use a new build directory when changing flags; make does not
track command-line flag changes. AddressSanitizer on macOS does not provide a
leak-check result; storage-release tests also check that the backend payload
pointer becomes null.

Tests do not read indeterminate values after allocation. Small integer-valued
real/complex fixtures are exactly representable, so these tests use
exact equality. Later numerical decomposition tests will use residual tolerances.
