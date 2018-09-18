include ../../$(PLAT).mk

SRC = $(wildcard $(FEBDIR)FEBio3/*.cpp)
OBJ = $(patsubst $(FEBDIR)FEBio3/%.cpp, %.o, $(SRC))
DEP = $(patsubst $(FEBDIR)FEBio3/%.cpp, %.d, $(SRC))

TARGET =  $(FEBDIR)build/bin/febio3.$(PLAT)

FELIBS =  $(FEBDIR)build/lib/libfecore_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebiolib_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebioplot_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebiomech_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebiomix_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebioxml_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libnumcore_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebioopt_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebiotest_$(PLAT).a
FELIBS += $(FEBDIR)build/lib/libfebiofluid_$(PLAT).a

FEBIOLIBS = -Wl,--start-group $(FELIBS) -Wl,--end-group

$(TARGET): $(OBJ) $(FELIBS)
ifeq ($(findstring lnx,$(PLAT)),lnx)
	$(CC) -o $(TARGET) $(DEF) $(FLG) -Wl,--export-dynamic $(INC) $(OBJ) $(FEBIOLIBS) $(LIBS) -ldl
else ifeq ($(findstring gcc,$(PLAT)),gcc)
	$(CC) -o $(TARGET) $(DEF) $(FLG) $(INC) $(OBJ) $(FEBIOLIBS) $(LIBS) -ldl
else ifeq ($(findstring sky,$(PLAT)),sky)
	$(CC) -o $(TARGET) $(DEF) $(FLG) $(INC) $(OBJ) $(FEBIOLIBS) $(LIBS) -ldl
else
	$(CC) -o $(TARGET) $(DEF) $(FLG) $(INC) $(OBJ) $(FELIBS) $(LIBS)
endif

%.o: $(FEBDIR)FEBio3/%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c $<

clean:
	$(RM) *.o *.d $(TARGET)

-include $(DEP)
