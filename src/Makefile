#
# Which compiler to use
#
#COMPILER = INTEL
COMPILER = GNU
#
#
ifeq "$(COMPILER)" "INTEL"
   MPIFORTRAN = mpiifort
   FORTRAN = ifort
   FLAGS = -O3
endif
ifeq "$(COMPILER)" "GNU"
   MPIFORTRAN = mpif90
   FORTRAN = gfortran
   FLAGS = -O3 -ffast-math
endif

# Directory for executables
BIN=../bin

ifeq ($(MAKECMDGOALS),devel)
  FLAGS = -Wall -fcheck=bounds
endif
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################

BINARIES = $(BIN)/compactlog 

MODULES = types.o

all : bindir $(BINARIES)
devel : bindir $(BINARIES)

#
# Link
#

$(BIN)/compactlog : Makefile $(MODULES) compactlog.o
	$(FORTRAN) -o $(BIN)/compactlog $(MODULES) compactlog.o $(FLAGS)

#
# Build objects
#

compactlog.o : Makefile compactlog.f90 
	$(FORTRAN) -c compactlog.f90 $(FLAGS)

#
# Build modules
#

modules : $(MODULES) 

types.o : types.f90
	$(FORTRAN) -c types.f90 $(FLAGS)

#
#  Binary directory
#

bindir : 
	@mkdir -p $(BIN)

#
# Clean everything
#

clean : 
	rm -f ./*.o ./*.mod

cleanall : 
	rm -f ./*.o ./*.mod $(BIN)/*





