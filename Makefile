# Compiler and flags
CC = g++
CFLAGS = -Wall -g -pthread
SINGLE = galign_benchmark
MULTI = multi_galign_benchmark
GPU = gpu_galign_benchmark
OBJDIR = build
INCLUDE_PATH = ./include
BUILD_PATH = ./build   # Compiled files will be placed here
LIBRARY_NAME = galign
LIBRARY_SRC = lib/cmd.cpp lib/fasta.cpp
LIBRARY_OBJ = $(OBJDIR)/cmd.o $(OBJDIR)/fasta.o
SHARED_OBJECT = $(OBJDIR)/lib$(LIBRARY_NAME).so

# Create build directories if they don't exist
$(shell mkdir -p $(OBJDIR) $(BUILD_PATH))

# Default target
all: $(SHARED_OBJECT) $(SINGLE) $(MULTI) $(GPU)

# Singlethreaded execution target
$(SINGLE): $(SINGLE).cpp $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -o $(SINGLE) $(SINGLE).cpp -I$(INCLUDE_PATH) -L$(BUILD_PATH) -lgalign 

# Multithreaded executable target
$(MULTI): $(MULTI).cpp $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -o $(MULTI) $(MULTI).cpp -I$(INCLUDE_PATH) -L$(BUILD_PATH) -lgalign 

# GPU executable target
$(GPU): $(GPU).cpp $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -o $(GPU) $(GPU).cpp -I$(INCLUDE_PATH) -L$(BUILD_PATH) -lgalign 

# Compile source files into object files in the build directory
$(OBJDIR)/%.o: lib/%.cpp
	$(CC) $(CFLAGS) -fPIC -c $< -I$(INCLUDE_PATH) -o $@

# Create shared library from object files
$(SHARED_OBJECT): $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -shared -o $@ $^

# Clean up generated files
clean:
	rm -rf $(BUILD_PATH)/lib$(LIBRARY_NAME).so $(SINGLE) $(MULTI) $(GPU) $(OBJDIR)

.PHONY: all clean
