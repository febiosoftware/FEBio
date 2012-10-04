LIB = ../lib/numcore_$(PLAT).a

$(LIB):
	$(CC) -c $(INC) $(DEF) $(FLG) *.cpp
	ar -cvr $(LIB) *.o

clean:
	rm *.o
	rm $(LIB)
