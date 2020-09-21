LIB = libjordanchevalley
CORE = $(LIB)/core
CPP = $(LIB)/cpp
PYTHON = $(LIB)/python

BUILD = build

CONDA_INCLUDE = conda-env/include
EIGEN_INCLUDE = conda-env/include/eigen3

CC = g++
FLAGS = -std=c++17 -Wall # -fsanitize=undefined -O3
INCLUDES = -I $(CONDA_INCLUDE) -I $(EIGEN_INCLUDE) -I $(CORE)


###################################################################
# COMPILE CPP BACKEND: make core
# COMPILE PYTHON FRONTEND: make py

# COMPILE CPP FRONTEND: make cpp
# RUN FRONTEND: make run
# COMPILE & RUN UNIT TESTS: make test

# CLEANUP: make clean
###################################################################

.PHONY: build cpp run py start 

build:
	mkdir -p $(BUILD)/$(CORE) $(BUILD)/$(CPP) $(BUILD)/$(PYTHON)

## Compile Object Files
# core -> build/core
$(BUILD)/$(CORE)/jordanchevalley.o: $(CORE)/jordanchevalley.cpp
	$(CC) $(FLAGS) $(INCLUDES) -c -o $@ $^ 

# cpp -> build/cpp
$(BUILD)/$(CPP)/io.o: $(CPP)/io.cpp
	$(CC) $(FLAGS) $(INCLUDES) -c -o $@ $^ 

$(BUILD)/$(CPP)/main.o: $(CPP)/main.cpp
	$(CC) $(FLAGS) $(INCLUDES) -c -o $@ $^ 


## Compile Linked Executables
# core shared library (needs compilation of all objects with -fPIC)
$(BUILD)/libjordanchevalley.so: $(CORE)/jordanchevalley.cpp
	$(CC) $(FLAGS) -mstackrealign -fPIC $(INCLUDES) --shared -o $@ $^

# cpp main.out executable
$(BUILD)/main.out: $(BUILD)/$(CPP)/main.o $(BUILD)/$(CPP)/io.o $(BUILD)/$(CORE)/jordanchevalley.o
	$(CC)  $(FLAGS) $(INCLUDES) -o $@ $^ 


## .PHONY targes for convenience
core: build $(BUILD)/$(CORE)/jordanchevalley.o
cpp: build $(BUILD)/main.out
run: cpp
	$(BUILD)/main.out


# build and run cpp test.out executable
test: build 
	$(CC) $(FLAGS) $(INCLUDES) -o $(BUILD)/test.out $(CPP)/unit_tests.cpp 
	$(BUILD)/test.out

# compile temporary python extension int ./build
py:
	python setup.py build_ext -b $(BUILD) -t $(BUILD) 

clean:
	rm -rf $(BUILD)