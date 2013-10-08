LIB = ../lib/febioxml_$(PLAT).a

$(LIB):
	$(CC) -c $(INC) $(DEF) $(FLG) *.cpp -fPIC
	ar -cvr $(LIB) *.o

clean:
	rm *.o
	rm $(LIB)
