
# source files
SRC := $(wildcard *.cpp)


OBJ = $(SRC:.cpp=.o)

OUT = simlib.a




#C++ compiler flags
CCFLAGS = -O2 -ow -std=c++17 -I$(shell  root-config --incdir)

#LD flags
LDFLAGS = -lgsl -lgslcblas $(shell root-config --libs)

#compiler
CCC = g++


.SUFFIXES: .cpp

default:  $(OUT)


.cpp.o:
	$(CCC) $(CCFLAGS) -c $< -o $@

$(OUT): $(OBJ)
	ar rcs $(OUT) $(OBJ)

clean:
	rm -f $(OBJ) $(OUT)
