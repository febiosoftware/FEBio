################################################################################
#                                                                              #
#                               Makefile for FEBio                             #
#                                                                              #
################################################################################

#  This can be whatever you want it to be, e.g., win for Windows, osx for Mac
PLAT = lnx

#  The c++ compiler, e.g. g++ for the gnu compiler
CC = g++

#  For Mac add: -DMAC, For Windows replace with: -DWIN32
DEF = -DLINUX

FLG = -static -O3

INC = -I..

# The default solver for FEBio is the Slyline solver, which is included with the
# FEBio source code.  Below are instructions for using the SuperLU and Pardiso solvers.
# SuperLU can be obtained from http://crd.lbl.gov/~xiaoye/SuperLU/.  It comes with its
# own BLAS routines, but it recommends using an optimized BLAS.  We are using the routines
# in the Intel Math Kernel Library (MKL).
# Pardiso can be obtained as a shared object library at http://www.pardiso-project.org
# We are using the Pardiso library contained in the MKL.

# Uncomment the following lines for the SuperLU solver
#DEF += -DSUPERLU
#  Put the path to the SuperLU directory here
#SUPERLU_PATH = 
#SUPERLU_LIB_PATH = $(SUPERLU_PATH)/lib
#  Make sure that the name of the library is correct
#SUPERLU_LIB = -L$(SUPERLU_LIB_PATH) -lsuperlu_3.1_$(PLAT)
#LIBS += $(SUPERLU_LIB)
#INC += -I$(SUPERLU_PATH)/SRC

# Both Pardiso and SuperLU use the Intel MKL.  Uncomment the following lines for the MKL
# Put the path to the MKL directory here
#MKL_PATH = 
#MKL_LIB = -L$(MKL_PATH) -lmkl_lapack -lmkl_em64t -liomp5 -lpthread
#LIBS += $(MKL_LIB)

# Uncomment the following line for the MKL Pardiso solver
#DEF += -DPARDISO


all:
	( cd FECore; $(MAKE) -f fecore.mk \
	  PLAT="$(PLAT)" CC="$(CC)" DEF="$(DEF)" FLG="$(FLG)" INC="$(INC)")
	( cd FEBio;  $(MAKE) -f febio.mk \
	  PLAT="$(PLAT)" CC="$(CC)" DEF="$(DEF)" FLG="$(FLG)" INC="$(INC)" LIBS="$(LIBS)")

clean:

	( cd FECore; $(MAKE) -f fecore.mk clean PLAT="$(PLAT)")
	( cd FEBio;  $(MAKE) -f febio.mk clean PLAT="$(PLAT)") 
