LIB = ../bin/febioopt_$(PLAT).so.1.0

FECORE = ../lib/fecore_$(PLAT).a

FEBIOXML = ../lib/febioxml_$(PLAT).a

FEBIO2LIBS = $(FEBIOXML) $(FECORE)

$(LIB):
	$(CC) -c $(INC) $(DEF) -O3 -Wall *.cpp -fPIC
	$(CC) -shared -Wl,-soname,febioopt_$(PLAT).so.1 -o ../bin/febioopt_$(PLAT).so.1.0 *.o $(FEBIO2LIBS)
	ln -sf ../bin/febioopt_$(PLAT).so.1.0 ../bin/febioopt_$(PLAT).so
	ln -sf ../bin/febioopt_$(PLAT).so.1.0 ../bin/febioopt_$(PLAT).so.1

clean:
	rm *.o
