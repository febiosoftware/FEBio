
TARGET = ../bin/febio2.$(PLAT)

FELIBS = ../lib/fecore_$(PLAT).a

FELIBS += ../lib/febiolib_$(PLAT).a

FELIBS += ../lib/febioplot_$(PLAT).a

FELIBS += ../lib/febiomech_$(PLAT).a

FELIBS += ../lib/febiomix_$(PLAT).a

FELIBS += ../lib/febioheat_$(PLAT).a

FELIBS += ../lib/febioxml_$(PLAT).a

FELIBS += ../lib/numcore_$(PLAT).a

FELIBS += ../lib/febioopt_$(PLAT).a

FELIBS += ../lib/liblevmar_$(PLAT).a

FEBIO2LIBS = -Wl,--start-group $(FELIBS) -Wl,--end-group

$(TARGET):
	$(CC) -o $(TARGET) $(DEF) *.cpp $(FLG) $(INC) $(FEBIO2LIBS) $(LIBS) -ldl

clean:
	rm -f $(TARGET)
