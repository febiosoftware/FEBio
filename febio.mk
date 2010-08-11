FEBIO = /home/sci/rawlins/FEBio
lnx: INCLUDE = $(FEBIO)/make.lnx
lnxdl: INCLUDE = $(FEBIO)/make.lnxdl
lnx32: INCLUDE = $(FEBIO)/make.lnx32
test: INCLUDE = $(FEBIO)/make.test
osx: INCLUDE = $(FEBIO)/make.osx
alt: INCLUDE = $(FEBIO)/make.alt
win: INCLUDE = $(FEBIO)/make.win
lnxclean: INCLUDE = $(FEBIO)/make.lnx
lnxdlclean: INCLUDE = $(FEBIO)/make.lnxdl
lnx32clean: INCLUDE = $(FEBIO)/make.lnx32
testclean: INCLUDE = $(FEBIO)/make.test
osxclean: INCLUDE = $(FEBIO)/make.osx
altclean: INCLUDE = $(FEBIO)/make.alt
winclean: INCLUDE = $(FEBIO)/make.win

export INCLUDE

lnx:
	./includes.bash
	./svnrev.bash
	( cd FECore; $(MAKE) )
	( cd FEBio; $(MAKE) )

lnxdl:
	./includes.bash
	./svnrev.bash
	( cd FECore; $(MAKE) )
	( cd FEBio; $(MAKE) )

lnx32:
	( cd FECore; $(MAKE) )
	( cd FEBio; $(MAKE) )

test:
	./includes.bash
	( cd FECore; $(MAKE) )
	( cd FEBio; $(MAKE) )

osx:
	( cd FECore; $(MAKE) )
	( cd FEBio; $(MAKE) )

alt:
	( cd FECore; $(MAKE) )
	( cd FEBio; $(MAKE) -f Makefile.alt )

win:
	( cd FECore; $(MAKE) -f Makefile.win )
	( cd FEBio; $(MAKE) -f Makefile.win )


lnxclean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 

lnxdlclean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 

lnx32clean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 

altclean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 

osxclean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 

winclean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 

testclean:

	( cd FECore; $(MAKE) clean)
	( cd FEBio;  $(MAKE) clean) 
