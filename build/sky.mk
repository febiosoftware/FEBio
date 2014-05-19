# Make include file for FEBio on Linux

CC = g++

DEF = -DLINUX

FLG = -O3 -fPIC

INC = -I$(FEBDIR)
