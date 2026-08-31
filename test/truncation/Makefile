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
truncation: clean main.o $(OBJECTS)
	$(CCCOM) $(CXX_FLAGS) -o truncation main.o $(OBJECTS) $(LIB_FLAGS)

clean:
	rm -f *.o

###############################################################
