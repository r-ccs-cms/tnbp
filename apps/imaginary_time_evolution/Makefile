# Supply the same CPU/CUDA configuration used by estimator.
CONFIGURATION ?= ../estimator/Configuration
include $(CONFIGURATION)
.DEFAULT_GOAL := imaginary_time_evolution
.PHONY: FORCE clean
imaginary_time_evolution: main.o
	$(CCCOM) $(CXX_FLAGS) -o $@ $^ $(LIB_FLAGS)
main.o: main.cc FORCE
	$(CCCOM) $(CXX_FLAGS) -c $< -o $@
FORCE:
clean:
	rm -f main.o
