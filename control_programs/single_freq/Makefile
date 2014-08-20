OS:=$(shell uname -s)

CCpp=g++
CFLAGS=-g
//LFLAGS=-L"${HOME}/uhd/host/build/lib"
libs=-lboost_program_options -lhdf5
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=single_freq

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CCpp) $(CFLAGS) -o $@ $(OBJECTS) $(libs)

.cpp.o:
	$(CCpp) $(CFLAGS) $< -o $@ -c

clean:
	rm -rf *.o client

