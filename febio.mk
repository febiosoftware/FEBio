FEBIO = $(CURDIR)
lnx64: INCLUDE = $(FEBIO)/make.lnx64
lnx64d: INCLUDE = $(FEBIO)/make.lnx64d
lnx32: INCLUDE = $(FEBIO)/make.lnx32
osx: INCLUDE = $(FEBIO)/make.osx
osxd: INCLUDE = $(FEBIO)/make.osxd
neohookeanpi: INCLUDE = $(FEBIO)/make.lnx64d
neohookeanpi_osx: INCLUDE = $(FEBIO)/make.osxd
lnx64clean: INCLUDE = $(FEBIO)/make.lnx64
lnx32clean: INCLUDE = $(FEBIO)/make.lnx32
osxclean: INCLUDE = $(FEBIO)/make.osx

export INCLUDE
export FEBIO

lnx64:
	( cd FEBioLib; $(MAKE) -f febiolib.mk )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk )
	( cd FEBioXML; $(MAKE) -f febioxml.mk )
	( cd NumCore; $(MAKE) -f numcore.mk )
	( cd FECore; $(MAKE) -f fecore.mk )
	( cd FEBioMix; $(MAKE) -f febiomix.mk )
	( cd FEBioMech; $(MAKE) -f febiomech.mk )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk )
	( cd FEBio2;  $(MAKE) -f febio2.mk )

lnx64d:
	( cd FEBioLib; $(MAKE) -f febiolib.mk )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk )
	( cd FEBioXML; $(MAKE) -f febioxml.mk )
	( cd NumCore; $(MAKE) -f numcore.mk )
	( cd FECore; $(MAKE) -f fecore.mk )
	( cd FEBioMix; $(MAKE) -f febiomix.mk )
	( cd FEBioMech; $(MAKE) -f febiomech.mk )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk )
	( cd FEBio2;  $(MAKE) -f febio2.mk )

lnx32:
	( cd FEBioLib; $(MAKE) -f febiolib.mk )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk )
	( cd FEBioXML; $(MAKE) -f febioxml.mk )
	( cd NumCore; $(MAKE) -f numcore.mk )
	( cd FECore; $(MAKE) -f fecore.mk )
	( cd FEBioMix; $(MAKE) -f febiomix.mk )
	( cd FEBioMech; $(MAKE) -f febiomech.mk )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk )
	( cd FEBio2;  $(MAKE) -f febio2.mk )

osx:
	( cd FEBioLib; $(MAKE) -f febiolib.mk )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk )
	( cd FEBioXML; $(MAKE) -f febioxml.mk )
	( cd NumCore; $(MAKE) -f numcore.mk )
	( cd FECore; $(MAKE) -f fecore.mk )
	( cd FEBioMix; $(MAKE) -f febiomix.mk )
	( cd FEBioMech; $(MAKE) -f febiomech.mk )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk )
	( cd FEBio2;  $(MAKE) -f febio2_osx.mk )

osxd:
	( cd FEBioLib; $(MAKE) -f febiolib.mk )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk )
	( cd FEBioXML; $(MAKE) -f febioxml.mk )
	( cd NumCore; $(MAKE) -f numcore.mk )
	( cd FECore; $(MAKE) -f fecore.mk )
	( cd FEBioMix; $(MAKE) -f febiomix.mk )
	( cd FEBioMech; $(MAKE) -f febiomech.mk )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk )
	( cd FEBio2;  $(MAKE) -f febio2_osx.mk )

neohookeanpi:
	( cd NeoHookeanPI; $(MAKE) -f neohookeanpi.mk )

neohookeanpi_osx:
	( cd NeoHookeanPI; $(MAKE) -f neohookeanpi_osx.mk )

lnx64clean:

	( cd FEBioLib; $(MAKE) -f febiolib.mk clean )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk clean )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk clean )
	( cd FEBioXML; $(MAKE) -f febioxml.mk clean )
	( cd NumCore; $(MAKE) -f numcore.mk clean )
	( cd FECore; $(MAKE) -f fecore.mk clean )
	( cd FEBioMix; $(MAKE) -f febiomix.mk clean )
	( cd FEBioMech; $(MAKE) -f febiomech.mk clean )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk clean )
	( cd FEBio2;  $(MAKE) -f febio2.mk clean )

lnx32clean:

	( cd FEBioLib; $(MAKE) -f febiolib.mk clean )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk clean )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk clean )
	( cd FEBioXML; $(MAKE) -f febioxml.mk clean )
	( cd NumCore; $(MAKE) -f numcore.mk clean )
	( cd FECore; $(MAKE) -f fecore.mk clean )
	( cd FEBioMix; $(MAKE) -f febiomix.mk clean )
	( cd FEBioMech; $(MAKE) -f febiomech.mk clean )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk clean )
	( cd FEBio2;  $(MAKE) -f febio2.mk clean )

osxclean:

	( cd FEBioLib; $(MAKE) -f febiolib.mk clean )
	( cd FEBioPlot; $(MAKE) -f febioplot.mk clean )
	( cd FEBioHeat; $(MAKE) -f febioheat.mk clean )
	( cd FEBioXML; $(MAKE) -f febioxml.mk clean )
	( cd NumCore; $(MAKE) -f numcore.mk clean )
	( cd FECore; $(MAKE) -f fecore.mk clean )
	( cd FEBioMix; $(MAKE) -f febiomix.mk clean )
	( cd FEBioMech; $(MAKE) -f febiomech.mk clean )
	( cd FEBioOpt; $(MAKE) -f febioopt.mk clean )
	( cd FEBio2;  $(MAKE) -f febio2.mk clean )
