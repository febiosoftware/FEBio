include $(FEBDIR)/make.$(PLAT)

SRC = $(wildcard ../*.cpp)
OBJ = $(patsubst ../%.cpp, %.o, $(SRC))
DEP = $(patsubst ../%.cpp, %.d, $(SRC))

SO = libneohookeanpi_$(PLAT).so
LIB = $(FEBDIR)/lib/$(SO)

FECORE = $(FEBDIR)/lib/libfecore_$(PLAT).a

FEBIOMECH = $(FEBDIR)/lib/libfebiomech_$(PLAT).a

FEBIOLIBS = $(FEBIOMECH) $(FECORE)

$(LIB): $(OBJ)
	$(CC) -shared -static-intel -no-intel-extensions -Wl,-soname,$(SO) -o $(LIB) $(OBJ) $(FEBIOLIBS)

%.o: ../%.cpp
	$(CC) $(INC) $(DEF) $(FLG) -MMD -c -o $@ $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
