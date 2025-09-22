include Configuration
TCI_INCLUDE_DIR=$(TCI_PATH)/include
TCL_INCLUDE_DIR=$(TCL_PATH)/include
HPTT_LIB=hptt
HPTT_INCLUDE_DIR=$(HPTT_PATH)/include
HPTT_LIB_DIR=$(HPTT_PATH)/lib
TNBP_INCLUDE_DIR=$(TNBP_PATH)/include
LIBFLAGS= -L$(HPTT_LIB_DIR) -l$(HPTT_LIB) $(SYSLIB)
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
	$(CCCOM) -c $(CCFLAGS) -I$(TCI_INCLUDE_DIR) -I$(TCL_INCLUDE_DIR) -I$(HPTT_INCLUDE_DIR) -I$(TNBP_INCLUDE_DIR) $<
###############################################################
kicked_ising_oned_exvalue: clean main.o $(OBJECTS)
	$(CCCOM) $(CCFLAGS) -o kicked_ising_oned_exvalue main.o $(OBJECTS) $(LIBFLAGS)

clean:
	rm -f *.o

###############################################################
