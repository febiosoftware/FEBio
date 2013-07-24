
TARGET = ../bin/febio2.$(PLAT)

FECORE = ../lib/fecore_$(PLAT).a

FEBIOLIB = ../lib/febiolib_$(PLAT).a

FEBIOPLOT = ../lib/febioplot_$(PLAT).a

FEBIOOPT = ../lib/febioopt_$(PLAT).a

FEBIOHEAT = ../lib/febioheat_$(PLAT).a

FEBIOXML = ../lib/febioxml_$(PLAT).a

NUMCORE = ../lib/numcore_$(PLAT).a

FEBIO2LIBS = -Wl,--start-group $(FEBIOLIB) $(FEBIOPLOT) $(FEBIOOPT) $(FEBIOHEAT) $(FEBIOXML) $(FECORE) $(NUMCORE) -Wl,--end-group

$(TARGET):
	$(CC) -o $(TARGET) $(DEF) *.cpp $(FLG) $(INC) $(FEBIO2LIBS) $(LIBS)

clean:
	rm $(TARGET)
