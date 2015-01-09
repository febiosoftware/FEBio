# Make include file for FEBio on Linux

include $(FEBDIR)build/lnx64.mk

CC = g++

FLG = -O3 -fPIC -fopenmp

INC = -I$(FEBDIR) -I$(FEBDIR)build/include
