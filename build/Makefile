FEBDIR = $(dir $(CURDIR))
lnx64:  PLAT = lnx64
lnx64d: PLAT = lnx64d
lnx64g: PLAT = lnx64g
lnx32:  PLAT = lnx32
gcc:    PLAT = gcc
osx:    PLAT = osx
osxd:   PLAT = osxd
sky:    PLAT = sky
lnx64clean:  PLAT = lnx64
lnx32clean:  PLAT = lnx32
osxclean:    PLAT = osx
lnx64dclean: PLAT = lnx64d
lnx64gclean: PLAT = lnx64g
gccclean:    PLAT = gcc
osxdclean:   PLAT = osxd
skyclean:    PLAT = sky

export PLAT
export FEBDIR

lnx64 lnx64d lnx64g lnx32 gcc osx osxd sky:
	( cd $(PLAT)/FEBioLib;  $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioPlot; $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioHeat; $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioXML;  $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/NumCore;   $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FECore;    $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioMix;  $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioMech; $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioOpt;  $(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBio2;    $(MAKE) -f ../../febio2.mk )

lnx64clean lnx32clean osxclean lnx64dclean lnx64gclean gccclean osxdclean skyclean:

	( cd $(PLAT)/FEBioLib;  $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioPlot; $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioHeat; $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioXML;  $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/NumCore;   $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FECore;    $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioMix;  $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioMech; $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioOpt;  $(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBio2;    $(MAKE) -f ../../febio2.mk clean )

.PHONY: lnx64 lnx32 osx lnx64d lnx64g gcc osxd sky