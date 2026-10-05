CONFIGURATION ?= Configuration
include $(CONFIGURATION)

.DEFAULT_GOAL := estimator
.PHONY: clean FORCE

estimator: main.o
	$(CCCOM) $(CXX_FLAGS) -o $@ $^ $(LIB_FLAGS)

# Rebuild when switching backend or command-line configuration. Keep this
# ordered before linking; a clean prerequisite would race with parallel builds.
main.o: main.cc FORCE
	$(CCCOM) $(CXX_FLAGS) -c $< -o $@

FORCE:

clean:
	rm -f main.o
