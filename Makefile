FEBDIR = $(CURDIR)
lnx64:  PLAT = lnx64
lnx64d: PLAT = lnx64d
lnx64g: PLAT = lnx64g
lnx32:  PLAT = lnx32
gcc:    PLAT = gcc
osx:    PLAT = osx
osxd:   PLAT = osxd
sky:  PLAT = sky
lnx64clean:  PLAT = lnx64
lnx32clean:  PLAT = lnx32
osxclean:    PLAT = osx
lnx64dclean: PLAT = lnx64d
lnx64gclean: PLAT = lnx64g
gccclean:    PLAT = gcc
osxdclean:   PLAT = osxd
skyclean:   PLAT = sky

export PLAT
export FEBDIR

lnx64 lnx64d lnx64g lnx32 gcc osx osxd sky:
	( cd build/$(PLAT)/FEBioLib;  $(MAKE) -f ../../Makelibs.mk FELIB="febiolib";  LIBDIR='FEBioLib')
	( cd build/$(PLAT)/FEBioPlot; $(MAKE) -f ../../Makelibs.mk FELIB="febioplot"; LIBDIR='FEBioPlot')
	( cd build/$(PLAT)/FEBioHeat; $(MAKE) -f ../../Makelibs.mk FELIB="febioheat"; LIBDIR='FEBioHeat')
	( cd build/$(PLAT)/FEBioXML;  $(MAKE) -f ../../Makelibs.mk FELIB="febioxml";  LIBDIR='FEBioXML')
	( cd build/$(PLAT)/NumCore;   $(MAKE) -f ../../Makelibs.mk FELIB="numcore";   LIBDIR='NumCore')
	( cd build/$(PLAT)/FECore;    $(MAKE) -f ../../Makelibs.mk FELIB="fecore";    LIBDIR='FECore')
	( cd build/$(PLAT)/FEBioMix;  $(MAKE) -f ../../Makelibs.mk FELIB="febiomix";  LIBDIR='FEBioMix')
	( cd build/$(PLAT)/FEBioMech; $(MAKE) -f ../../Makelibs.mk FELIB="febiomech"; LIBDIR='FEBioMech')
	( cd build/$(PLAT)/FEBioOpt;  $(MAKE) -f ../../Makelibs.mk FELIB="febioopt";  LIBDIR='FEBioOpt')
	( cd build/$(PLAT)/FEBio2;    $(MAKE) -f ../../febio2.mk )

lnx64clean lnx32clean osxclean lnx64dclean lnx64gclean gccclean osxdclean skyclean:

	( cd build/$(PLAT)/FEBioLib;  $(MAKE) -f ../../Makelibs.mk FELIB="febiolib" clean )
	( cd build/$(PLAT)/FEBioPlot; $(MAKE) -f ../../Makelibs.mk FELIB="febioplot" clean )
	( cd build/$(PLAT)/FEBioHeat; $(MAKE) -f ../../Makelibs.mk FELIB="febioheat" clean )
	( cd build/$(PLAT)/FEBioXML;  $(MAKE) -f ../../Makelibs.mk FELIB="febioxml" clean )
	( cd build/$(PLAT)/NumCore;   $(MAKE) -f ../../Makelibs.mk FELIB="numcore" clean )
	( cd build/$(PLAT)/FECore;    $(MAKE) -f ../../Makelibs.mk FELIB="fecore" clean )
	( cd build/$(PLAT)/FEBioMix;  $(MAKE) -f ../../Makelibs.mk FELIB="febiomix" clean )
	( cd build/$(PLAT)/FEBioMech; $(MAKE) -f ../../Makelibs.mk FELIB="febiomech" clean )
	( cd build/$(PLAT)/FEBioOpt;  $(MAKE) -f ../../Makelibs.mk FELIB="febioopt" clean )
	( cd build/$(PLAT)/FEBio2;    $(MAKE) -f ../../febio2.mk clean )

.PHONY: lnx64 lnx32 osx lnx64d lnx64g gcc osxd sky