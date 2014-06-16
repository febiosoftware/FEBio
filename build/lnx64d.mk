# Make include file for FEBio on Linux 64 bit

CC = icpc

# Remove -DHAVE_LEVMAR and $(LEV_LIB) from LIBS if not linking with the Lourakis levmar routine.
DEF = -DLINUX -DPARDISO -DHAVE_LEVMAR -DHAVE_ZLIB

FLG = -O3 -fPIC -openmp -static-intel -no-intel-extensions

#Intel Compiler
INTELROOT = $(subst /mkl,,$(MKLROOT))/compiler
INTEL_INC = $(INTELROOT)/include
INTEL_LIB = $(INTELROOT)/lib/intel64

# Pardiso solver
MKL_PATH = $(MKLROOT)/lib/intel64
MKL_LIB = -Wl,--start-group $(MKL_PATH)/libmkl_intel_lp64.a
MKL_LIB += $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group
MKL_LIB += $(INTEL_LIB)/libiomp5.a -pthread -lz

#Levmar library
LEV_LIB = -llevmar_$(PLAT)

LIBS = -L$(FEBDIR)/build/lib $(LEV_LIB) $(MKL_LIB)

INC = -I$(INTEL_INC) -I$(FEBDIR) -I$(FEBDIR)build/include
