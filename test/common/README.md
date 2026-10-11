# Shared test fixtures

These helpers have no standalone test registration or product implementation.

- `graph_fixture.h`: deterministic real/complex data, independent dense
  open-boundary enumeration, and shared numerical/metadata checks for generic
  graph relocation and ordered BMPS sweeps. The sweep operations are supplied
  as callbacks; the fixture does not choose a BMPS policy.
- `dense_fixture.h`: deterministic tensor-entry formula shared by the MCMC
  reference calculation and the BMPS upper-region independence test.

Test cases and executable entry points remain in their owning suite directories.
