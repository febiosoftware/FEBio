LIB = ../bin/febioopt_$(PLAT).so.1.0


$(LIB):
	$(CC) -c $(INC) $(DEF) $(FLG) *.cpp -Wall -fPIC
	$(CC) -shared -Wl,-soname,febioopt_$(PLAT).so.1 -o ../bin/febioopt_$(PLAT).so.1.0 *.o
	ln -sf ../bin/febioopt_$(PLAT).so.1.0 ../bin/febioopt_$(PLAT).so
	ln -sf ../bin/febioopt_$(PLAT).so.1.0 ../bin/febioopt_$(PLAT).so.1

clean:
	rm *.o
