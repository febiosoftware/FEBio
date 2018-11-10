# Make include file for FEBio on Mac

CC = icpc

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
DEF = -DPARDISO -DMKL_ISS -DHAVE_LEVMAR -DHAVE_GSL -DSVN

# Comment this out during development
DEF += -DNDEBUG

FLG = -O3 -fopenmp -fPIC -static-intel -no-intel-extensions -std=c++11

# Pardiso solver
INTELROOT = $(subst /mkl,,$(MKLROOT))
INTEL_INC = $(INTELROOT)/compiler/include
INTEL_LIB = $(INTELROOT)/compiler/lib/
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_intel_thread.a $(MKL_PATH)libmkl_core.a \
	$(INTEL_LIB)libiomp5.a -pthread -lz

#Levmar library
LEV_LIB = -llevmar

GSL_LIB = /usr/local/lib/libgsl.a

LIBS = -L$(FEBDIR)build/lib $(LEV_LIB) $(GSL_LIB) $(MKL_LIB)

INC = -I$(INTEL_INC) -I$(FEBDIR) -I$(FEBDIR)build/include
