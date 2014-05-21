# Make include file for FEBio on Linux

include $(FEBDIR)build/make.lnx64d

FLG := $(FLG:O3=g) # Note that we had to use := so that FLG is not recursive

