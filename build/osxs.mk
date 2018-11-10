# Make include file for FEBio sequential on Mac

CC = icpc

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
DEF = -DPARDISO -DHAVE_LEVMAR -DHAVE_ZLIB -DSVN

FLG = -O3 -fPIC -static-intel -no-intel-extensions -std=c++11

# Pardiso solver
INTELROOT = $(subst /mkl,,$(MKLROOT))
INTEL_INC = $(INTELROOT)/compiler/include
INTEL_LIB = $(INTELROOT)/compiler/lib/
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_sequential.a $(MKL_PATH)libmkl_core.a \
	-pthread -lz -lm

#Levmar library
LEV_LIB = -llevmar_$(PLAT)

LIBS = -L$(FEBDIR)build/lib $(LEV_LIB) $(MKL_LIB)

INC = -I$(INTEL_INC) -I$(FEBDIR) -I$(FEBDIR)build/include
