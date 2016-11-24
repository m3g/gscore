#
# Which compiler to use
#
FORTRAN = gfortran

FLAGS = -O3 --fast-math
ifeq ($(MAKECMDGOALS),devel)
  FLAGS = -Wall -fcheck=bounds
endif
 
# Directory for executables
BIN=../bin

###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################

BINARIES = $(BIN)/compactlog \
           $(BIN)/xcompactlog \
           $(BIN)/gscore \
           $(BIN)/qmatrix \
           $(BIN)/xgscore \
           $(BIN)/gcorrelation \
           $(BIN)/xgcorrelation \
           $(BIN)/xcompare \
           $(BIN)/dgscore \
           $(BIN)/evalconstraints \
           $(BIN)/checkconstraints

MODULES = types.o \
          compactlog_data.o \
          file_operations.o

all : $(BIN) $(BINARIES)
devel : $(BIN) $(BINARIES)

#
#  Binaries directory
#

$(BIN) : 
	@mkdir -p $(BIN)


#
# Link
#

$(BIN)/compactlog : Makefile $(MODULES) compactlog.o \
                                        model_index.o \
                                        sort.o \
                                        read_compactlog.o \
                                        title.o \
                                        progress_bar.o
	$(FORTRAN) -o $(BIN)/compactlog $(MODULES) \
                    compactlog.o \
                    model_index.o \
                    sort.o \
                    read_compactlog.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/xcompactlog : Makefile $(MODULES) xcompactlog.o \
                                         model_index.o \
                                         sort.o \
                                         read_compactlog.o \
                                         title.o \
                                         progress_bar.o
	$(FORTRAN) -o $(BIN)/xcompactlog $(MODULES) \
                    xcompactlog.o \
                    model_index.o \
                    sort.o \
                    read_compactlog.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/gscore : Makefile $(MODULES) gscore.o \
                                    sort.o \
                                    read_compactlog.o \
                                    title.o \
                                    progress_bar.o
	$(FORTRAN) -o $(BIN)/gscore $(MODULES) \
                    gscore.o \
                    sort.o \
                    read_compactlog.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/qmatrix : Makefile $(MODULES) qmatrix.o \
                                    title.o \
                                    sort.o \
                                    progress_bar.o
	$(FORTRAN) -o $(BIN)/qmatrix $(MODULES) \
                    qmatrix.o \
                    sort.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/xgscore : Makefile $(MODULES) xgscore.o \
                                     sort.o \
                                     read_compactlog.o \
                                     title.o \
                                     progress_bar.o
	$(FORTRAN) -o $(BIN)/xgscore $(MODULES) \
                    xgscore.o \
                    sort.o \
                    read_compactlog.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/gcorrelation : Makefile $(MODULES) gcorrelation.o \
                                          model_index.o \
                                          sort.o \
                                          read_compactlog.o \
                                          title.o \
                                          progress_bar.o
	$(FORTRAN) -o $(BIN)/gcorrelation $(MODULES) \
                    gcorrelation.o \
                    model_index.o \
                    sort.o \
                    read_compactlog.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/xgcorrelation : Makefile $(MODULES) xgcorrelation.o \
                                           model_index.o \
                                           sort.o \
                                           read_compactlog.o \
                                           title.o \
                                           progress_bar.o
	$(FORTRAN) -o $(BIN)/xgcorrelation $(MODULES) \
                    xgcorrelation.o \
                    model_index.o \
                    sort.o \
                    read_compactlog.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/xcompare : Makefile $(MODULES) xcompare.o \
                                      model_index.o \
                                      sort.o \
                                      title.o \
                                      progress_bar.o
	$(FORTRAN) -o $(BIN)/xcompare $(MODULES) \
                    xcompare.o \
                    model_index.o \
                    sort.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/dgscore : Makefile $(MODULES) dgscore.o \
                                    model_index.o \
                                    sort.o \
                                    title.o \
                                    progress_bar.o
	$(FORTRAN) -o $(BIN)/dgscore $(MODULES) \
                    dgscore.o \
                    model_index.o \
                    sort.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/evalconstraints : Makefile $(MODULES) evalconstraints.o \
                                    model_index.o \
                                    sort.o \
                                    title.o \
                                    progress_bar.o
	$(FORTRAN) -o $(BIN)/evalconstraints $(MODULES) \
                    evalconstraints.o \
                    model_index.o \
                    sort.o \
                    title.o \
                    progress_bar.o \
                    $(FLAGS)

$(BIN)/checkconstraints : Makefile $(MODULES) checkconstraints.o \
                                    model_index.o \
                                    sort.o \
                                    title.o \
                                    progress_bar.o
	$(FORTRAN) -o $(BIN)/checkconstraints $(MODULES) \
                    checkconstraints.o \
                    model_index.o \
                    sort.o \
                    title.o \
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

qmatrix.o : Makefile qmatrix.f90 
	$(FORTRAN) -c qmatrix.f90 $(FLAGS)

xgscore.o : Makefile xgscore.f90 
	$(FORTRAN) -c xgscore.f90 $(FLAGS)

gcorrelation.o : Makefile gcorrelation.f90 
	$(FORTRAN) -c gcorrelation.f90 $(FLAGS)

xgcorrelation.o : Makefile xgcorrelation.f90 
	$(FORTRAN) -c xgcorrelation.f90 $(FLAGS)

xcompare.o : Makefile xcompare.f90 
	$(FORTRAN) -c xcompare.f90 $(FLAGS)

dgscore.o : Makefile dgscore.f90 
	$(FORTRAN) -c dgscore.f90 $(FLAGS)

evalconstraints.o : Makefile evalconstraints.f90 
	$(FORTRAN) -c evalconstraints.f90 $(FLAGS)

checkconstraints.o : Makefile checkconstraints.f90 
	$(FORTRAN) -c checkconstraints.f90 $(FLAGS)

#
# Subroutine and function objects
#

model_index.o : Makefile model_index.f90 
	$(FORTRAN) -c model_index.f90 $(FLAGS)

progress_bar.o : Makefile progress_bar.f90 
	$(FORTRAN) -c progress_bar.f90 $(FLAGS)

title.o : Makefile title.f90 
	$(FORTRAN) -c title.f90 $(FLAGS)

sort.o : Makefile sort.f90 
	$(FORTRAN) -c sort.f90 $(FLAGS)

read_compactlog.o : Makefile read_compactlog.f90 
	$(FORTRAN) -c read_compactlog.f90 $(FLAGS)

#
# Build modules
#

modules : $(MODULES) 

compactlog_data.o : compactlog_data.f90 
	$(FORTRAN) -c compactlog_data.f90 $(FLAGS)

types.o : types.f90
	$(FORTRAN) -c types.f90 $(FLAGS)

file_operations.o : file_operations.f90
	$(FORTRAN) -c file_operations.f90 $(FLAGS)

#
# Clean everything
#

clean : 
	rm -f ./*.o ./*.mod

cleanall : 
	rm -f ./*.o ./*.mod $(BINARIES)





