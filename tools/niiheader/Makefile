SHELL = /bin/sh
# specify flags and libraries needed for your compiler
CC=mpicxx
CCFLAGS=-g -O3 -std=c++11
LDFLAGS=
OBJ=niiheader.o nifti1_io.o znzlib.o
EXE=niiheader

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp nifti1.h nifti1_io.h znzlib.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -rf $(OBJ) $(EXE) main.o
