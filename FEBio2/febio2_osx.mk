include $(FEBDIR)/make.$(PLAT)

SRC = $(wildcard ../*.cpp)
OBJ = $(patsubst ../%.cpp, %.o, $(SRC))
DEP = $(patsubst ../%.cpp, %.d, $(SRC))

TARGET =  $(FEBDIR)/bin/febio2.$(PLAT)

FELIBS =  $(FEBDIR)/lib/libfecore_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebiolib_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebioplot_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebiomech_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebiomix_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebioheat_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebioxml_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libnumcore_$(PLAT).a
FELIBS += $(FEBDIR)/lib/libfebioopt_$(PLAT).a
FELIBS += $(FEBDIR)/lib/liblevmar_$(PLAT).a

$(TARGET): $(OBJ)
	$(CC) -o $(TARGET) $(DEF) $(FLG) $(INC) $(OBJ) $(FELIBS) $(LIBS)

%.o: ../%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(TARGET)

-include $(DEP)
