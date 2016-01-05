# Makefile

INSTALL_DIR       := ./

INCLUDES          += -I. -I/usr/include/mpich
LIBS              += -L. -lmpi -lgsl -lgslcblas

CC              := g++

################################################################################

################################################################################
# Rules for Makefile
################################################################################
all: DDpar DRpar RRpar

DDpar : DDpar.c
	$(CXX) DDpar.c -o DDpar $(INCLUDES) $(LIBS)

DRpar : DRpar.c
	$(CXX) DRpar.c -o DRpar $(INCLUDES) $(LIBS)

RRpar : RRpar.c
	$(CXX) RRpar.c -o RRpar $(INCLUDES) $(LIBS)

clean:
	rm -f $(INSTALL_DIR)/DDpar
	rm -f $(INSTALL_DIR)/DRpar
	rm -f $(INSTALL_DIR)/RRpar

################################################################################
