# Make include file for FEBio on Mac
# This makefile assumes that llvm was installed with homebrew in its default directories (brew install llvm)

LLVM_PATH = /usr/local/opt/llvm

CC = $(LLVM_PATH)/bin/clang++

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
# Remove -DHAVE_GSL and $(GSL_LIB) from LIBS if not linking with the GNU scientific library.
DEF = -DPARDISO -DMKL_ISS -DHAVE_LEVMAR -DHAVE_GSL -DSVN -DHYPRE -DUSE_MPI -DHAS_MMG -DNDEBUG

FLG = -O3 -fopenmp -fPIC -std=c++11

# Pardiso solver
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_intel_thread.a $(MKL_PATH)libmkl_core.a

#Levmar library
LEV_PATH = /usr/local
LEV_LIB = -llevmar

#GSL library
GSL_PATH = /usr/local/opt/gsl
GSL_LIB = -lgsl

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

LIBS = -L$(LLVM_PATH)/lib -L$(FEBDIR)build/lib -L$(HYPRE_PATH)/lib $(HYPRE_LIB) -L$(LEV_PATH)/lib $(LEV_LIB) -L$(GSL_PATH)/lib $(GSL_LIB) $(MKL_LIB) -L$(OMP_PATH)/lib $(OMP_LIB) -L$(MPI_PATH)/lib $(MPI_LIB) -L$(MMG_PATH)/lib $(MMG_LIB)

INC = -I$(LLVM_PATH)/include -I$(FEBDIR) -I$(FEBDIR)build/include -I$(HYPRE_PATH)/include -I$(OMP_PATH)/include -I$(MPI_PATH)/include -I$(MMG_PATH)/include -I$(GSL_PATH)/include -I$(LEV_PATH)/include
