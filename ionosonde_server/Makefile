OS:=$(shell uname -s)

CCpp=g++ -std=c++0x
CC=gcc
CFLAGS=-c -g3
common_libs=-lm -lrt -luhd -lfftw3 -lboost_system -lboost_thread -lboost_filesystem
#LFLAGS=-L"${HOME}/uhd/host/build/lib"
INCLUDES=-I"${HOME}/uhd/host/include" -I"${HOME}/uhd/host/include/uhd/types" \
		 -I"${HOME}/uhd/host/include/uhd/usrp" -I"." -I"/usr/include/boost"
CPP_SOURCES=main.cpp capture_spectrum.cpp \
			lp_filter.cpp doppler_process.cpp \
			matched_filter.cpp transceive.cpp utils.cpp
CPP_OBJECTS=$(CPP_SOURCES:.cpp=.o)
#EXECUTABLE=../build/ionosonde_server
EXECUTABLE=ionosonde_server

all: $(EXECUTABLE)

$(EXECUTABLE): $(CPP_OBJECTS)
	$(CCpp) -o $@ $(CPP_OBJECTS) $(linux_libs) $(common_libs)

.cpp.o:
	$(CCpp) $(CFLAGS) $< -o $@ -c $(INCLUDES) $(common_libs)

clean:
	rm -rf *.o sounder

