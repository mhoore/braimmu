SHELL = /bin/sh
# specify flags and libraries needed for your compiler
CC=mpicxx
CCFLAGS=-g -O3 -std=c++11
LDFLAGS=
OBJ=main.o scenario_connectome.o scenario_geometry.o input.o init.o comm.o output.o region.o nifti1_io.o znzlib.o strategy.o
EXE=braimmu.exe

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp pointers.h scenario_connectome.h scenario_geometry.h input.h init.h comm.h output.h region.h nifti1.h nifti1_io.h znzlib.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -rf $(OBJ) $(EXE) main.o
