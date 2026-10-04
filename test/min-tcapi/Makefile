# Supply platform-specific compiler, OpenMP, MPI and BLAS flags externally.
CXXFLAGS ?= -std=c++17 -O2
ROOT := ../..
CPPFLAGS += -I$(ROOT)/external/min-tci/include -I$(ROOT)/external/tensor-ng-dev/Source/gqten/include -I$(ROOT)/external/tensor-ng-dev/Source/gqten/ext/hptt/include
BUILD_DIR ?= build
GROUPS := type_system read_only_queries construction_and_destruction miscellaneous tensor_manipulation tensor_linear_algebra io diagnostics
BINS := $(addprefix $(BUILD_DIR)/,$(GROUPS))
HEADERS := $(shell find $(ROOT)/external/min-tci/include/tcapi . -name '*.h')
.PHONY: all check
all: $(BINS)
$(BUILD_DIR):
	mkdir -p $@
$(BUILD_DIR)/type_system: type_system/main.cc type_system/other.cc $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) type_system/main.cc type_system/other.cc $(LDFLAGS) $(LDLIBS) -o $@
$(BUILD_DIR)/%: %/main.cc $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
check: all
	@set -e; for test in $(BINS); do echo "$$test"; "$$test"; done

# Optional TNBP integration target; run with mpirun -np 1 and -np 2.
# Kept outside check because it requires an MPI runtime.
$(BUILD_DIR)/tnbp: CPPFLAGS += -I$(ROOT)/include
