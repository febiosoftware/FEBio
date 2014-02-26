FEBDIR = $(CURDIR)
lnx64: PLAT = lnx64
lnx64d: PLAT = lnx64d
lnx32: PLAT = lnx32
osx: PLAT = osx
osxd: PLAT = osxd
neohookeanpi: PLAT = lnx64d
neohookeanpi_osx: PLAT = osxd
lnx64clean: PLAT = lnx64
lnx32clean: PLAT = lnx32
osxclean: PLAT = osx

export PLAT
export FEBDIR

lnx64 lnx64d lnx32:
	( cd FEBioLib/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiolib" )
	( cd FEBioPlot/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioplot" )
	( cd FEBioHeat/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioheat" )
	( cd FEBioXML/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioxml" )
	( cd NumCore/$(PLAT);   $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="numcore" )
	( cd FECore/$(PLAT);    $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="fecore" )
	( cd FEBioMix/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiomix" )
	( cd FEBioMech/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiomech" )
	( cd FEBioOpt/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioopt" )
	( cd FEBio2/$(PLAT);    $(MAKE) -f ../febio2.mk )

osx osxd:
	( cd FEBioLib/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiolib" )
	( cd FEBioPlot/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioplot" )
	( cd FEBioHeat/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioheat" )
	( cd FEBioXML/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioxml" )
	( cd NumCore/$(PLAT);   $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="numcore" )
	( cd FECore/$(PLAT);    $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="fecore" )
	( cd FEBioMix/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiomix" )
	( cd FEBioMech/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiomech" )
	( cd FEBioOpt/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioopt" )
	( cd FEBio2/$(PLAT);    $(MAKE) -f ../febio2_osx.mk )

neohookeanpi:
	( cd NeoHookeanPI/$(PLAT); $(MAKE) -f ../neohookeanpi.mk )

neohookeanpi_osx:
	( cd NeoHookeanPI/$(PLAT); $(MAKE) -f ../neohookeanpi_osx.mk )

lnx64clean lnx32clean osxclean:

	( cd FEBioLib/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiolib" clean )
	( cd FEBioPlot/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioplot" clean )
	( cd FEBioHeat/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioheat" clean )
	( cd FEBioXML/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioxml" clean )
	( cd NumCore/$(PLAT);   $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="numcore" clean )
	( cd FECore/$(PLAT);    $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="fecore" clean )
	( cd FEBioMix/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiomix" clean )
	( cd FEBioMech/$(PLAT); $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febiomech" clean )
	( cd FEBioOpt/$(PLAT);  $(MAKE) -f $(FEBDIR)/Makelibs.mk FELIB="febioopt" clean )
	( cd FEBio2/$(PLAT);    $(MAKE) -f ../febio2.mk clean )
