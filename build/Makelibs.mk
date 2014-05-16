include $(FEBDIR)/make.$(PLAT)

SRC = $(wildcard ../*.cpp)
OBJ = $(patsubst ../%.cpp, %.o, $(SRC))
DEP = $(patsubst ../%.cpp, %.d, $(SRC))

LIB = $(FEBDIR)/lib/lib$(FELIB)_$(PLAT).a

$(LIB): $(OBJ)
	ar -cvr $(LIB) $(OBJ)

%.o: ../%.cpp
	$(CC) -MMD -c $(INC) $(DEF) $(FLG) $<

clean:
	$(RM) *.o *.d $(LIB)

-include $(DEP)
