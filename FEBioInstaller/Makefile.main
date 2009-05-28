all:
	( cd FECore; $(MAKE) )
	( cd FEBio;  $(MAKE) )

linux:
	cp make.lnx make.inc
	( cd FECore; $(MAKE) clean; $(MAKE) )
	( cd FEBio; $(MAKE) clean;  $(MAKE) )

mac:
	cp make.osx make.inc
	( cd FECore; $(MAKE) clean; $(MAKE) )
	( cd FEBio; $(MAKE) clean;  $(MAKE) )

altix:
	cp make.alt make.inc
	( cd FECore; $(MAKE) clean; $(MAKE) )
	( cd FEBio; $(MAKE) clean;  $(MAKE) -f Makefile.alt )


clean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 
