OS:=$(shell uname -s)

CCpp=g++
LFLAGS = -L"${HOME}/uhd/host/build/lib"
common_libs=-lm -lrt -luhd -lfftw3
linux_libs=-lboost_thread -lboost_program_options
INCLUDES=-I"${HOME}/uhd/host/include" -I"${HOME}/uhd/host/include/uhd/types" -I"${HOME}/uhd/host/include/uhd/usrp" -I"."
CPP_SOURCES=main.cpp recv_clr_freq.cpp recv_to_file.cpp tx_worker.cpp
CPP_OBJECTS=$(CPP_SOURCES:.cpp=.o)
EXECUTABLE=sounder

all: $(EXECUTABLE)

$(EXECUTABLE): $(CPP_OBJECTS)
	$(CCpp) -o $@ $(CPP_OBJECTS) $(linux_libs) $(common_libs)

.cpp.o:
	$(CCpp) $(CFLAGS) $< -o $@ -c $(LFLAGS) $(INCLUDES)

clean:
	rm -rf *.o sounder

