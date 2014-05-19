# Make include file for FEBio on Mac
# Debug version (doesn't include -DNDEBUG)

CC = icpc

DEF = -DPARDISO -DHAVE_LEVMAR

FLG = -O3 -fopenmp -fPIC -static-intel -no-intel-extensions

# Pardiso solver
INTELROOT = $(subst /mkl,,$(MKLROOT))
INTEL_INC = $(INTELROOT)/compiler/include
INTEL_LIB = $(INTELROOT)/compiler/lib/
MKL_PATH = $(MKLROOT)/lib/
MKL_LIB = $(MKL_PATH)libmkl_intel_lp64.a $(MKL_PATH)libmkl_intel_thread.a $(MKL_PATH)libmkl_core.a \
	$(MKL_PATH)libmkl_pgi_thread.a $(INTEL_LIB)libiomp5.a

# KeyGen
KEYGEN = -lkeygen_$(PLAT)

#Levmar library
LEV_LIB = -llevmar_$(PLAT)

LIBS = -L$(FEBDIR)build/lib $(LEV_LIB) $(MKL_LIB)

INC = -I$(INTEL_INC) -I$(FEBDIR) -I$(FEBDIR)include -I/usr/include/c++/4.2.1
