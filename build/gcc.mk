# Make include file for FEBio on Linux

include $(FEBDIR)build/lnx64d.mk

CC = g++

FLG = -O3 -fPIC

INC = -I$(FEBDIR) -I$(FEBDIR)build/include
