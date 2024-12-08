# Compiler and flags
CC = g++
NVCC = nvcc
CFLAGS = -Wall -g -pthread
CUDAFLAGS = -g

SINGLE = galign_benchmark
MULTI = multi_galign_benchmark
GPU = gpu_galign_benchmark
OBJDIR = build
INCLUDE_PATH = ./include
BUILD_PATH = ./build   # Shared object and object files will be placed here
LIBRARY_NAME = galign
LIBRARY_SRC = lib/cmd.cpp lib/fasta.cpp lib/globals.cpp
LIBRARY_OBJ = $(OBJDIR)/cmd.o $(OBJDIR)/fasta.o $(OBJDIR)/globals.o
SHARED_OBJECT = $(OBJDIR)/lib$(LIBRARY_NAME).so

# Create build directory if it doesn't exist
$(shell mkdir -p $(OBJDIR) $(BUILD_PATH))

# Default target
all: $(SHARED_OBJECT) $(SINGLE) $(MULTI) $(GPU)

# Library target
lib: $(SHARED_OBJECT)

# Library and Singlethreaded target
single: lib $(SINGLE)

# Library and Multithreaded target
multi: lib $(MULTI)

# Library and GPU target
gpu: lib $(GPU)

# Compile source files into object files in the build directory
$(OBJDIR)/%.o: lib/%.cpp
	$(CC) $(CFLAGS) -fPIC -c $< -I$(INCLUDE_PATH) -o $@

# Create shared library from object files
$(SHARED_OBJECT): $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -shared -o $@ $^

# Singlethreaded execution target
$(SINGLE): $(SINGLE).cpp $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -o $(SINGLE) $(SINGLE).cpp -I$(INCLUDE_PATH) -L$(BUILD_PATH) -lgalign 

# Multithreaded execution target
$(MULTI): $(MULTI).cpp $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -o $(MULTI) $(MULTI).cpp -I$(INCLUDE_PATH) -L$(BUILD_PATH) -lgalign 

# GPU execution target
$(GPU): $(GPU).cu $(LIBRARY_OBJ)
	$(NVCC) $(CUDAFLAGS) -o $(GPU) $(GPU).cu -I$(INCLUDE_PATH) -L$(BUILD_PATH) -lgalign 

# Clean up generated files
clean:
	rm -rf $(BUILD_PATH)/lib$(LIBRARY_NAME).so $(SINGLE) $(MULTI) $(GPU) $(OBJDIR)

.PHONY: all lib single multi gpu clean
