################################################################################
#                                                                              #
#                               Makefile for FEBio                             #
#                                                                              #
################################################################################

#  This can be whatever you want it to be, e.g., win for Windows, osx for Mac
PLAT = lnx

#  The c++ compiler, e.g. g++ for the gnu compiler
CC = g++

#  For Mac add: -DMAC. For Windows replace with: -DWIN32
DEF = -DLINUX

FLG = -O3

INC = -I..

# The default solver for FEBio is the Slyline solver, which is included with the
# FEBio source code.  Below are instructions for using the SuperLU and Pardiso solvers.
# SuperLU can be obtained from http://crd.lbl.gov/~xiaoye/SuperLU/.  It comes with its
# own BLAS routines, but it recommends using an optimized BLAS.
# Pardiso can be obtained as a shared object library at http://www.pardiso-project.org
# Another version of Pardiso is included in the Intel MKL.  These two version are no
# longer compatible.

# Uncomment the following lines for the SuperLU solver (version 4.0)
# SuperLU comes with a version of the BLAS library, but we recommend using a BLAS optimized
# for your platform.  We are using the BLAS in the MKL.
#DEF += -DSUPERLU
#  Put the path to the SuperLU directory here
#SUPERLU_PATH = 
#SUPERLU_LIB_PATH = $(SUPERLU_PATH)/lib
#  Make sure that the name of the library is correct
#SUPERLU_LIB = -L$(SUPERLU_LIB_PATH) -lsuperlu_4.0_$(PLAT)
#LIBS += $(SUPERLU_LIB)
#INC += -I$(SUPERLU_PATH)/SRC

# Uncomment the following lines for the shared object library version of the Pardiso solver.
# Note that you can use either the MKL version or the .so version, but not both.
# This version on linux also requires a fortran library and an optimized BLAS. We are using
# the Intel Fortran library and the BLAS in the MKL.
# You may need to change this for your implementation.
#DEF += -DPARDISODL
# Put the path to the Pardiso library here
#PARDISO_PATH = 
# Put the path to the Fortran library here
#FORT_PATH = 
# We are using the BLAS from the MKL, so uncomment the MKL lines below (but not the line with
# -DPARDISO) or add directions to your version of the BLAS.
#PARDISO_LIB = -L$(PARDISO_PATH) -lpardiso400_INTEL101_IA64 -L$(FORT_PATH) -lifcore
#LIBS += $(PARDISO_LIB)

# Uncomment the following lines for the MKL Pardiso solver
# Linking Pardiso from the MKL depends on your machine architecture and the version of the MKL.
# This example is for the MKL 10.2.3 on a 64bit linux machine.  See the MKL documentation
# for your architecture and MKL version.
#DEF += -DPARDISO
# Put the path to the MKL library directory here
#MKL_PATH = 
#MKL_LIB = -L$(MKL_PATH) -lmkl_intel_lp64
#MKL_LIB += -Wl,--start-group -lmkl_intel_thread -lmkl_core -Wl,--end-group
#MKL_LIB += -liomp5 -lpthread
#LIBS += $(MKL_LIB)


all:
	( cd FECore; $(MAKE) -f fecore.mk \
	  PLAT="$(PLAT)" CC="$(CC)" DEF="$(DEF)" FLG="$(FLG)" INC="$(INC)")
	( cd FEBio;  $(MAKE) -f febio.mk \
	  PLAT="$(PLAT)" CC="$(CC)" DEF="$(DEF)" FLG="$(FLG)" INC="$(INC)" LIBS="$(LIBS)")

clean:

	( cd FECore; $(MAKE) -f fecore.mk clean PLAT="$(PLAT)")
	( cd FEBio;  $(MAKE) -f febio.mk clean PLAT="$(PLAT)") 
