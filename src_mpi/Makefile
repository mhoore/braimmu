SHELL = /bin/sh
# specify flags and libraries needed for your compiler
CC=mpicxx
CCFLAGS=-g -O3 -std=c++14
LDFLAGS=
OBJ=main.o memory.o brain.o input.o init.o comm.o run.o output.o region.o nifti1_io.o znzlib.o NiftiView.o
EXE=braimmu.exe

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp pointers.h memory.h brain.h input.h init.h comm.h run.h output.h region.h nifti1.h nifti1_io.h znzlib.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -rf $(OBJ) $(EXE) main.o
