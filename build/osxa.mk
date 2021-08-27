# Make include file for FEBio on Mac using Apple compiler

LLVM_PATH = /usr

CC = $(LLVM_PATH)/bin/clang++

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
DEF = -DPARDISO -DMKL_ISS -DHAVE_LEVMAR -DSVN -DHYPRE -DUSE_MPI -DHAS_MMG -DNDEBUG

FLG = -Ofast -Xpreprocessor -fopenmp -fPIC -std=c++17

# Pardiso solver
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_intel_thread.a $(MKL_PATH)libmkl_core.a

#Levmar library
LEV_PATH = /usr/local
LEV_LIB = -llevmar

#HYPRE library
HYPRE_LIB = -lHYPRE
HYPRE_PATH = /usr/local/hypre-master/src/hypre

#OMP library
OMP_PATH = /usr/local/opt/libomp
OMP_LIB = -lomp

#MPI library
MPI_PATH = /usr/local/opt/open-mpi
MPI_LIB = -lmpi

#MMG library
MMG_PATH = /usr/local/mmg-master/build
MMG_LIB = -lmmg3d

LIBS = -L$(LLVM_PATH)/lib -L$(FEBDIR)build/lib -L$(HYPRE_PATH)/lib $(HYPRE_LIB) -L$(LEV_PATH)/lib $(LEV_LIB) $(MKL_LIB) -L$(OMP_PATH)/lib $(OMP_LIB) -L$(MPI_PATH)/lib $(MPI_LIB) -L$(MMG_PATH)/lib $(MMG_LIB)

INC = -I$(MKLROOT)/include -I$(LLVM_PATH)/include -I$(FEBDIR) -I$(FEBDIR)build/include -I$(HYPRE_PATH)/include -I$(OMP_PATH)/include -I$(MPI_PATH)/include -I$(MMG_PATH)/include -I$(LEV_PATH)/include
