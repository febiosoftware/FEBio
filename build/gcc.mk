# Make include file for FEBio on Linux

include $(FEBDIR)build/lnx32.mk

CC = g++

FLG = -O3 -fPIC -fopenmp -std=c++11

INC = -I$(FEBDIR) -I$(FEBDIR)build/include
