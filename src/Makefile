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
   FLAGS = -O3
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

BINARIES = $(BIN)/compactlog \
           $(BIN)/xcompactlog \
           $(BIN)/gscore \
           $(BIN)/gcorrelation

MODULES = types.o \
          file_operations.o

all : bindir $(BINARIES)
devel : bindir $(BINARIES)

#
# Link
#

$(BIN)/compactlog : Makefile $(MODULES) compactlog.o \
                                        model_index.o \
                                        sort.o \
                                        progress_bar.o
	$(FORTRAN) -o $(BIN)/compactlog $(MODULES) \
                    compactlog.o \
                    model_index.o \
                    sort.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/xcompactlog : Makefile $(MODULES) xcompactlog.o \
                                        model_index.o \
                                        sort.o \
                                        progress_bar.o
	$(FORTRAN) -o $(BIN)/xcompactlog $(MODULES) \
                    xcompactlog.o \
                    model_index.o \
                    sort.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/gscore : Makefile $(MODULES) gscore.o \
                                    sort.o \
                                    progress_bar.o
	$(FORTRAN) -o $(BIN)/gscore $(MODULES) \
                    gscore.o \
                    sort.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/gcorrelation : Makefile $(MODULES) gcorrelation.o \
                                        model_index.o \
                                        sort.o \
                                        progress_bar.o
	$(FORTRAN) -o $(BIN)/gcorrelation $(MODULES) \
                    gcorrelation.o \
                    model_index.o \
                    sort.o \
                    progress_bar.o \
                    $(FLAGS)


#
# Main program objects
#

compactlog.o : Makefile compactlog.f90 
	$(FORTRAN) -c compactlog.f90 $(FLAGS)

xcompactlog.o : Makefile xcompactlog.f90 
	$(FORTRAN) -c xcompactlog.f90 $(FLAGS)

gscore.o : Makefile gscore.f90 
	$(FORTRAN) -c gscore.f90 $(FLAGS)

gcorrelation.o : Makefile gcorrelation.f90 
	$(FORTRAN) -c gcorrelation.f90 $(FLAGS)

#
# Subroutine and function objects
#

model_index.o : Makefile model_index.f90 
	$(FORTRAN) -c model_index.f90 $(FLAGS)

progress_bar.o : Makefile progress_bar.f90 
	$(FORTRAN) -c progress_bar.f90 $(FLAGS)

sort.o : Makefile sort.f90 
	$(FORTRAN) -c sort.f90 $(FLAGS)

#
# Build modules
#

modules : $(MODULES) 

types.o : types.f90
	$(FORTRAN) -c types.f90 $(FLAGS)

file_operations.o : file_operations.f90
	$(FORTRAN) -c file_operations.f90 $(FLAGS)

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





