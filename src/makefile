# Compiler
CXX = mpiCC

# Compiler flags
CXXFLAGS = -g -O3 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -lm -fopenmp -Wall -Werror 

# Source files
SRCS = main.cpp circuit.cpp common.cpp ini.cpp INIReader.cpp simulator.cpp IO_Runner.cpp DIO_Runner.cpp MEM_Runner.cpp MPI_Runner.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = Quokka

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

# Rule to compile the source files
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
