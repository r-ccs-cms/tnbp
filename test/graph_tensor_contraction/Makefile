include Configuration
MAKEFILES= Makefile Configuration
# header file
HEADER=
# source
SOURCES= main.cc
#objects
OBJECTS=
# compilation
.SUFFIXES:
.SUFFIXES: .h .cc .o
.cc.o:  $*.cc
	$(CCCOM) -c $(CXX_FLAGS) $<
###############################################################
gtnc: clean main.o $(OBJECTS)
	$(CCCOM) $(CXX_FLAGS) -o graph_tensor_contraction main.o $(OBJECTS) $(LIB_FLAGS)

clean:
	rm -f *.o

###############################################################
