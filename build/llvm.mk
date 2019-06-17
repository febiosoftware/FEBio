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
LEV_LIB = -llevmar

#GSL library
GSL_LIB = -lgsl

#HYPRE library
HYPRE_LIB = -lHYPRE
HYPRE_PATH = /usr/local/hypre-2.11.2/src/hypre

#OMP library
OMP_PATH = /usr/local/opt/libomp

#MPI library
MPI_LIB = -lmpi

#MMG library
MMG_LIB = -lmmg3d

LIBS = -L$(LLVM_PATH)/lib -L$(FEBDIR)build/lib -L$(HYPRE_PATH)/lib $(HYPRE_LIB) $(LEV_LIB) $(GSL_LIB) $(MKL_LIB) $(MPI_LIB) $(MMG_LIB)

INC = -I$(LLVM_PATH)/include -I$(FEBDIR) -I$(FEBDIR)build/include -I$(HYPRE_PATH)/include -I$(OMP_LIB)/include
