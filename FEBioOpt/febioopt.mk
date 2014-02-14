LIB = ../lib/febioopt_$(PLAT).a

$(LIB):
	$(CC) -c $(INC) $(DEF) $(FLG) -DHAVE_LEVMAR *.cpp
	ar -cvr $(LIB) *.o

clean:
	rm -f *.o
	rm -f $(LIB)
