# Compiler and flags
CC = g++
CFLAGS = -Wall -g
LIBRARY_NAME = galignlib
LIBRARY_SRC = lib/cmd.cpp lib/parser.cpp
LIBRARY_OBJ = lib/cmd.o lib/parser.o
SINGLE = galign_benchmark
MULTI = galign_multi
GPU = galign_gpu
LIBRARY_PATH = /lib
INCLUDE_PATH = /include

# Targets
all: $(SINGLE) # $(MULTI) $(GPU)

$(SINGLE): galign_benchmark.cpp $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -o $(SINGLE) galign_benchmark.cpp -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -l$(LIBRARY_NAME)

# Not yet implemented
# $(MULTI): galign_multi.cpp $(LIBRARY_OBJ)
# 	$(CC) $(CFLAGS) -o $(EXEC) galign_multi.cpp -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -l$(LIBRARY_NAME)

# $(GPU): galign_gpu.cpp $(LIBRARY_OBJ)
# 	$(CC) $(CFLAGS) -o $(EXEC) galign_gpu.cpp -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -l$(LIBRARY_NAME)

$(LIBRARY_OBJ): $(LIBRARY_SRC)
	$(CC) $(CFLAGS) -fPIC -c $(LIBRARY_SRC) -I$(INCLUDE_PATH)

$(LIBRARY_NAME): $(LIBRARY_OBJ)
	$(CC) $(CFLAGS) -shared -o $(LIBRARY_PATH)/lib$(LIBRARY_NAME).so $(LIBRARY_OBJ)

clean:
	rm -f $(LIBRARY_PATH)/lib$(LIBRARY_NAME).so $(SINGLE) $(LIBRARY_OBJ)
