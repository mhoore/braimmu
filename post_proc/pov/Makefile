SHELL = /bin/sh
# specify flags and libraries needed for your compiler
CC=g++
CCFLAGS=-g -O3 -std=c++11
LDFLAGS=
OBJ=nii2pov.o nifti1_io.o znzlib.o
EXE=n2p.exe

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.cpp nii2pov.h nifti1.h nifti1_io.h znzlib.h
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -rf $(OBJ) $(EXE) main.o
