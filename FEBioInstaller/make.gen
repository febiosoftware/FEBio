# Make include file for FEBio on Linux

HOME_PATH = /home/sci/rawlins/Projects

PLAT = storm

CC = g++

DEF = -DLINUX -DSUPERLU -DPARDISO

FLG = -O3

# Pardiso solver
FORTRAN_LIB = -L/usr/lib64/gcc/x86_64-suse-linux/4.2.1/ -lgfortran
ACML = /opt/acml4.2.0/gfortran64/lib
ACML_LIB = $(ACML)/libacml.a $(ACML)/libacml_mv.a -lgomp
PARDISO_LIB = -L$(HOME_PATH)/Pardiso/ -lpardiso_GNU42_AMD64_INT64_P

# SuperLU solver
SUPERLU_INC = $(HOME_PATH)/SuperLU_3.1/SRC/
SUPERLU_LIB_PATH = $(HOME_PATH)/SuperLU_3.1/lib/
SUPERLU_LIB = -L$(SUPERLU_LIB_PATH) -lsuperlu_3.1_$(PLAT) -lblas

LIBS = $(SUPERLU_LIB) $(PARDISO_LIB) $(ACML_LIB) $(FORTRAN_LIB)

INC = -I$(SUPERLU_INC) -I..
