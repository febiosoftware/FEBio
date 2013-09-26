FEBIO = /home/sci/rawlins/FEBio2
lnx64: INCLUDE = $(FEBIO)/make.lnx64
lnx64clean: INCLUDE = $(FEBIO)/make.lnx64

export INCLUDE

lnx64:
	( cd FEBioPlot; $(MAKE) )
	( cd FEBioLib; $(MAKE) )
	( cd FEBioOpt; $(MAKE) )
	( cd FEBioXML; $(MAKE) )
	( cd NumCore; $(MAKE) )
	( cd FECore; $(MAKE) )
	( cd FEBio2; $(MAKE) )



lnx64clean:

	( cd FEBioPlot; $(MAKE) clean)
	( cd FEBioLib; $(MAKE) clean)
	( cd FEBioOpt; $(MAKE) clean)
	( cd FEBioXML; $(MAKE) clean)
	( cd NumCore; $(MAKE) clean)
	( cd FECore; $(MAKE) clean)
	( cd FEBio2; $(MAKE) clean) 

