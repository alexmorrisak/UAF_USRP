OS:=$(shell uname -s)

CCpp=g++ -std=c++0x
CC=gcc
CFLAGS=-c -g
LFLAGS=-L"${HOME}/uhd/host/build/lib"
common_libs=-lm -lrt -luhd -lfftw3
linux_libs=-lboost_thread -lboost_program_options
INCLUDES=-I"${HOME}/uhd/host/include" -I"${HOME}/uhd/host/include/uhd/types" \
		 -I"${HOME}/uhd/host/include/uhd/usrp" -I"."
CPP_SOURCES=main.cpp capture_spectrum.cpp \
			lp_filter.cpp doppler_process.cpp \
			matched_filter.cpp transceive.cpp utils.cpp
C_SOURCES=
CPP_OBJECTS=$(CPP_SOURCES:.cpp=.o)
C_OBJECTS=$(C_SOURCES:.c=.o)
EXECUTABLE=../build/ionosonde_server

all: $(EXECUTABLE)

$(EXECUTABLE): $(CPP_OBJECTS) $(C_OBJECTS)
	$(CCpp) -o $@ $(CPP_OBJECTS) $(C_OBJECTS) $(linux_libs) $(common_libs) $(LFLAGS)

.cpp.o:
	$(CCpp) $(CFLAGS) $< -o $@ -c $(LFLAGS) $(INCLUDES)

.c.o:
	$(CC) $(CFLAGS) $< -o $@ -c $(LFLAGS) $(INCLUDES)

clean:
	rm -rf *.o sounder

