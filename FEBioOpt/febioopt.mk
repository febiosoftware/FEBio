SO = libfebioopt_$(PLAT).so
LIB = ../lib/$(SO)

FECORE = ../lib/fecore_$(PLAT).a

FEBIOXML = ../lib/febioxml_$(PLAT).a

FEBIO2LIBS = $(FEBIOXML) $(FECORE)

$(LIB):
	$(CC) -c $(INC) $(DEF) -O3 -Wall *.cpp -fPIC
	$(CC) -shared -Wl,-soname,$(SO).1 -o $(LIB).1.0 *.o $(FEBIO2LIBS)
	ln -sf $(LIB).1.0 $(LIB)
	ln -sf $(LIB).1.0 $(LIB).1

clean:
	rm -f *.o
