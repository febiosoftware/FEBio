# ************************************  FEBio Makefile  *************************************
#
# The file Mkdir.bash needs to be run to create the directories of the configurations
# you want to build.  For example, execute "./Mkdir.bash lnx64" to create the directories
# for lnx64.  Then run "make lnx64" to build that configuration.
# Several configurations are included in this directory.  All configurations except "sky"
# assume that the Intel MKL is installed and that the environment variable MKLROOT is defined
# (running # the Intel script "compilervars" will set the necessary environment variables).
# sky: Compiles using g++ and uses the default skyline linear solver.  You will need to
# edit febio.xml in the bin directory (replace "pardiso" with "skyline") to use this solver.

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