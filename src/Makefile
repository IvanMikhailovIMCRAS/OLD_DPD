MF=      Makefile
FC=      gfortran
FFLAGS=  -O3 -fopenmp
LFLAGS=  $(FFLAGS)

# name of EXEcutable file:
EXE=     ../dpd.exe

# list of SouRCe files:
SRC= \
         CommonParam.f90  \
         ErrorList.f90 \
         BoxMuller.f90 \
         Lib.f90 \
         Output.f90 \
         ForceEnergy.f90 \
         Engine.f90 \
         DPD.f90

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)
.f90.o:
	$(FC) $(FFLAGS) -c $<
all:	$(EXE)
$(EXE):	$(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)
$(OBJ):	$(MF)
tar:
	tar cvf $(EXE).tar $(MF) $(SRC)
clean:
	rm -f $(OBJ) $(EXE) *.mod core
