
TARGET = ../bin/febio.$(PLAT)

FECORE = ../lib/fecore_$(PLAT).a

$(TARGET):
	$(CC) -o $(TARGET) $(DEF) $(INC) *.cpp $(FLG) $(FECORE) $(LIBS)

clean:
	rm $(TARGET)
