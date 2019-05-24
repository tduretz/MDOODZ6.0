# MDOODZ 6.0 makefile

# Compiler
CC = gcc

#---------------------------------------------------#
# C flags
CFLAGS = -std=c99 -D _UMFPACK_ -I ./

ifeq ($(OPT),yes)	
	CFLAGS += -O3 -ftree-vectorize -funroll-loops
else
	CFLAGS += -g -Wall  -O0 -fno-inline -fno-omit-frame-pointer
endif

ifeq ($(OMP),yes)
	CFLAGS += -fopenmp -D _OMP_
else
	CFLAGS += -Wno-unknown-pragmas 
endif

ifeq ($(VG),yes)
	CFLAGS += -D _VG_ -Wno-format-zero-length
endif

CFLAGS += -Wno-unused-variable -Wno-unused-but-set-variable -Wno-comment
#---------------------------------------------------#
# Libraries
LIBS = -lz -lhdf5

# Related to SuiteSparse
LIBS += -llapack -lcxsparse -lumfpack -lcblas -lamd -lcholmod -lcolamd -lbtf -lsuitesparseconfig -lblas

# Link to openmp
ifeq ($(OMP),yes)
    LIBS += -lgomp
endif

#---------------------------------------------------#
# Rules
all: Doodzi_$(MODEL)

#FILES = Main_MDOODZ.o ./$(MODEL_PATH)/set_$(MODEL).o StokesRoutines.o StokesAssemblyCoupled.o StokesAssemblyDecoupled.o AdvectionRoutines.o RheologyDensity.o OutputHDF5.o SparseTools.o ThermalRoutines.o ThermalSolver.o ParticleRoutines.o free_surf.o FlowLaws.o MemoryAllocFree.o InputOutput.o MiscFunctions.o GridRoutines.o

#FILES = StokesAssemblyCoupled_BEN.o Routines_BEN.o free_surf_BEN.o Main_MDOODZ_BEN_dec.o ./$(MODEL_PATH)/set_$(MODEL).o StokesRoutines.o StokesAssemblyCoupled.o StokesAssemblyDecoupled.o AdvectionRoutines.o RheologyDensity.o OutputHDF5.o SparseTools.o ThermalRoutines.o ThermalSolver.o ParticleRoutines.o free_surf.o FlowLaws.o MemoryAllocFree.o InputOutput.o MiscFunctions.o GridRoutines.o

FILES = Main_MDOODZ_GREAT.o Routines_BEN.o ./$(MODEL_PATH)/set_$(MODEL).o StokesRoutines.o StokesAssemblyCoupled.o StokesAssemblyDecoupled.o AdvectionRoutines.o RheologyDensity.o OutputHDF5.o SparseTools.o ThermalRoutines.o ThermalSolver.o ParticleRoutines.o free_surf.o FlowLaws.o MemoryAllocFree.o InputOutput.o MiscFunctions.o GridRoutines.o BC.o User.o

Doodzi_$(MODEL): ${FILES}
	$(CC) ${FILES} -o Doodzi_$(MODEL) ${LIBS}

clean:
	rm -rf *o ./$(MODEL_PATH)/*.o Doodzi_$(MODEL)
