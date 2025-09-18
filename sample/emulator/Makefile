include Configuration
TCI_INCLUDE_DIR=$(TCI_PATH)/include
TCL_INCLUDE_DIR=$(TCL_PATH)/include
AEL_INCLUDE_DIR=$(AEL_PATH)/include
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
	$(CCCOM) -c $(CCFLAGS) -I$(TCI_INCLUDE_DIR) -I$(TCL_INCLUDE_DIR) -I$(AEL_INCLUDE_DIR) -I$(TNBP_INCLUDE_DIR) $<
###############################################################
test: clean main.o $(OBJECTS)
	$(CCCOM) $(CCFLAGS) -o test main.o $(OBJECTS) $(LIBFLAGS)

clean:
	rm -f *.o

###############################################################
