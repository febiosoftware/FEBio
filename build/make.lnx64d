# Make include file for FEBio on Linux

CC = icpc

# Add DPRINTHB for matrix output, DNAGLIB for NAG optimization
# Add  DHAVE_LEVMAR for Levmar optimization
# Compile with DNDEBUG for release version
#DEF = -DLINUX -DSUPERLU -DPARDISO -DHAVE_LEVMAR -DSVN -DFEBIO_LICENSE
DEF = -DLINUX -DPARDISO -DHAVE_LEVMAR

FLG = -O3 -fPIC -openmp -static-intel -no-intel-extensions

#Intel Compiler
INTELROOT = $(subst /mkl,,$(MKLROOT))/compiler
INTEL_INC = $(INTELROOT)/include
INTEL_LIB = $(INTELROOT)/lib/intel64

# Pardiso solver
MKL_PATH = $(MKLROOT)/lib/intel64
MKL_LIB = -Wl,--start-group $(MKL_PATH)/libmkl_intel_lp64.a
MKL_LIB += $(MKL_PATH)/libmkl_intel_thread.a $(MKL_PATH)/libmkl_core.a -Wl,--end-group
MKL_LIB += $(INTEL_LIB)/libiomp5.a -pthread

# KeyGen
KEYGEN = -lkeygen_$(PLAT)

#NAG library
NAG_INC=/usr/sci/apps/mrl/NAG/cll6a09dhl/include
NAG_LIB=-L/usr/sci/apps/mrl/NAG/cll6a09dhl/lib -lnagc_nag -lpthread

#Levmar library
LEV_LIB = -llevmar_$(PLAT)

LIBS = -L$(FEBDIR)/build/lib $(LEV_LIB) $(MKL_LIB)

INC = -I$(INTEL_INC) -I$(FEBDIR) -I$(FEBDIR)build/include
