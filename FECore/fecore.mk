LIB = ../lib/fecore_$(PLAT).a

$(LIB):
	$(CC) -c $(INC) $(DEF) $(FLG) *.cpp
	ar -cvq $(LIB) *.o

clean:
	rm *.o
	rm $(LIB)
