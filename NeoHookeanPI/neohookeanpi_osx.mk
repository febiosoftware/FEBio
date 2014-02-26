include $(FEBDIR)/make.$(PLAT)

SRC = $(wildcard ../*.cpp)
OBJ = $(patsubst ../%.cpp, %.o, $(SRC))
DEP = $(patsubst ../%.cpp, %.d, $(SRC))

LIB = $(FEBDIR)/lib/libneohookeanpi_$(PLAT).dylib

FECORE = $(FEBDIR)/lib/libfecore_$(PLAT).a

FEBIOMECH = $(FEBDIR)/lib/libfebiomech_$(PLAT).a

FEBIOLIBS = $(FEBIOMECH) $(FECORE)

$(LIB): $(OBJ)
	$(CC) -dynamiclib $(FLG) $(INC) $(DEF) -o $(LIB) $(OBJ) $(FEBIOLIBS)

%.o: ../%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
