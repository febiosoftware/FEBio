# Make include file for FEBio on Mac

CC = clang++

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
DEF = -DPARDISO -DHAVE_LEVMAR -DSVN

FLG = -O3 -fopenmp -fPIC -std=c++11 -stdlib=libstdc++

# openMP
OMP_PATH = /usr/local/
OMP_LIB = $(OMP_PATH)/lib/libiomp5.dylib

# Pardiso solver
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_intel_thread.a $(MKL_PATH)libmkl_core.a

#Levmar library
LEV_LIB = -llevmar

LIBS = -L$(FEBDIR)build/lib $(LEV_LIB) $(MKL_LIB) $(OMP_LIB)

INC = -I$(FEBDIR) -I$(FEBDIR)build/include -I$(OMP)/include
