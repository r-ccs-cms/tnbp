# TNBP integration checks

This optional MPI test checks the migrated TNBP call sites: real/complex
truncation with a known spectrum on one or two ranks, diagonal messenger shapes,
TPS stream save/load, explicit TensorProductState copying, and SquareRoot.

Build the `build/tnbp` target of the parent Makefile with the same compiler and
library flags as the adapter tests (use an MPI compiler). Run from the parent:

```
mpirun -np 1 build/tnbp
mpirun -np 2 build/tnbp
```

It is intentionally outside `make check`, which does not require MPI launching.
The existing `test/truncation` and `test/square_root_and_inverse` provide the
larger numerical regression fixtures.
