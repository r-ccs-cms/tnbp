include Configuration
TNBP_INCLUDE_DIR=$(TNBP_PATH)/include
LIBFLAGS= $(SYSLIB)
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
	$(CCCOM) -c $(CCFLAGS) -I$(TNBP_INCLUDE_DIR) $<
###############################################################
parser: clean main.o $(OBJECTS)
	$(CCCOM) $(CCFLAGS) -o parser main.o $(OBJECTS) $(LIBFLAGS)

clean:
	rm -f *.o

###############################################################
