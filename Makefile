# --- Variable definitions ---
CC = icc
CXX = icc
#CC = g++

#CFLAGS = -Wall -std=c++20 -I./include 
CPPFLAGS = -qopenmp -Wall -I"${MKLROOT}/include" -fPIC -DMKL_ILP64 -m64 
CXXFLAGS = -std=c++17 
CPPFLAGS += -fnon-call-exceptions

#LDFLAGS = -lm -lrt -fopenmp -fexceptions
LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lm -ldl

OBJFILES = fdtd.o helpers.o
HFILES = include/helpers.h include/parameters.h

TARGET = test.out

# --- Rule 1 ---
all: clean $(TARGET)

$(TARGET): $(OBJFILES) $(HFILES)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
