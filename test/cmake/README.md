# Common test build configuration

`Backend.cmake` creates the `tnbp_test_backend` interface target and selects
CPU gqten or tcapi-cuda dependencies. It also places executables in the build
root. Both the preferred `test/CMakeLists.txt` entry and the legacy
`test/backend_integration/CMakeLists.txt` entry include this file.

Numerical test targets and CTest registrations belong to the individual suite
directories; see [the test index](../README.md).
