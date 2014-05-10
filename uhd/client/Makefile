OS:=$(shell uname -s)

CCpp=g++
//CFLAGS=
//LFLAGS=-L"${HOME}/uhd/host/build/lib"
libs=-lboost_program_options -lhdf5
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=client

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CCpp) -o $@ $(OBJECTS) $(libs)

.cpp.o:
	$(CCpp) $< -o $@ -c

clean:
	rm -rf *.o client

