LIB = ../lib/numcore_$(PLAT).a

$(LIB):
	$(CC) -c $(INC) $(DEF) $(FLG) *.cpp -fPIC
	ar -cvr $(LIB) *.o

clean:
	rm -f *.o
	rm -f $(LIB)
