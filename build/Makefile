# ************************************  FEBio Makefile  *************************************
#
# The file Mkdir.bash needs to be run to create the directories of the configurations
# you want to build.  For example, execute "./Mkdir.bash lnx64" to create the directories
# for lnx64.  Then run "make lnx64" to build that configuration.
# Several configurations are included in this directory.  All configurations except "sky"
# assume that the Intel MKL is installed and that the environment variable MKLROOT is defined
# (running # the Intel script "compilervars" will set the necessary environment variables).
# They also assume that you are linking with the Lourakis levmar routine, which can be
# downloaded at http://users.ics.forth.gr/~lourakis/levmar/.
# The main configuration files are lnx64d, osxd, lnx32 and sky.  The other are derived from
# those.
# sky: Compiles using g++ and uses the default skyline linear solver.  You will need to
# edit febio.xml in the bin directory (replace "pardiso" with "skyline") to use this solver.

FEBDIR = $(dir $(CURDIR))
lnx64:  PLAT = lnx64
lnx64d: PLAT = lnx64d
lnx64g: PLAT = lnx64g
lnx64s: PLAT = lnx64s
lnx64gs: PLAT = lnx64gs
lnx32:  PLAT = lnx32
gcc:    PLAT = gcc
gcc64:    PLAT = gcc64
gcc64g:    PLAT = gcc64g
osx:    PLAT = osx
osxd:   PLAT = osxd
osxs:   PLAT = osxs
clang: PLAT = clang
llvm: PLAT = llvm
sky:    PLAT = sky
lnx64clean:  PLAT = lnx64
lnx32clean:  PLAT = lnx32
osxclean:    PLAT = osx
lnx64dclean: PLAT = lnx64d
lnx64gclean: PLAT = lnx64g
lnx64sclean: PLAT = lnx64s
lnx64gsclean: PLAT = lnx64gs
gccclean:    PLAT = gcc
gcc64clean:    PLAT = gcc64
gcc64gclean:    PLAT = gcc64g
osxdclean:   PLAT = osxd
osxsclean:   PLAT = osxs
clangclean:  PLAT = clang
llvmclean:  PLAT = llvm
skyclean:    PLAT = sky

export PLAT
export FEBDIR

lnx64 lnx64d lnx64g lnx64s lnx64gs lnx32 gcc gcc64 gcc64g osx osxd osxs clang llvm sky:
	( cd $(PLAT)/FEBioLib;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioPlot;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioXML;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/NumCore;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FECore;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioMix;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioMech;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioOpt;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioTest;		$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBioFluid;	$(MAKE) -f ../../Makelibs.mk )
	( cd $(PLAT)/FEBio3;		$(MAKE) -f ../../febio3.mk )

lnx64clean lnx32clean osxclean lnx64dclean lnx64gclean lnx64sclean lnx64gsclean gccclean gcc64clean gcc64gclean osxdclean osxsclean clangclean llvmclean skyclean:
	( cd $(PLAT)/FEBioLib;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioPlot;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioXML;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/NumCore;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FECore;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioMix;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioMech;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioOpt;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioTest;		$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBioFluid;	$(MAKE) -f ../../Makelibs.mk clean )
	( cd $(PLAT)/FEBio3;		$(MAKE) -f ../../febio3.mk clean )


.PHONY: lnx64 lnx32 osx lnx64d lnx64g lnx64s lnx64gs gcc gcc64 gcc64g osxd osxs clang llvm sky
