# Make include file for FEBio on Linux

include $(FEBDIR)build/lnx64d.mk

FLG := $(FLG:O3=g) # Note that we had to use := so that FLG is not recursive
