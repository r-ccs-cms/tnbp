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
ptnos-bp: clean main.o $(OBJECTS)
	$(CCCOM) $(CXX_FLAGS) -o ptnos-bp main.o $(OBJECTS) $(LIB_FLAGS)

clean:
	rm -f *.o

###############################################################
