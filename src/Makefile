##===============================================##
##  Code to Generate and Analyze Rough Surfaces  ##
##                                               ##
##  V.A. Yastrebov, 2013-2015                    ##
##  Centre des Materiaux, CNRS, MINES ParisTech  ##
##===============================================##

CC = g++ 
## You'll need for define a path to FFTW3
INCDIR = -I$(FFTPATH)/include
INCLIB = -L$(FFTPATH)/lib
OPTIM    =  -O2 # -g 
WARNFLAG = -Wall
FULLFLAG = $(OPTIM) $(WARNFLAG) $(INCDIR) $(INCLIB) #-std=gnu++0x -DUSING_IOHELPER #-DUSING_QVIEW
CFLAGS = -lfftw3 -lstdc++
FILE_C = $(wildcard *.c)
FILE_GENERATOR = generate_random_surface.c 
EXE = SURFACE_GENERATOR

all: generate_random_surface 

generate_random_surface: generate_random_surface.c normal.c
	$(CC) $(FULLFLAG) $^ -o $(EXE) $(CFLAGS); mv $(EXE) ../bin
#$(CC) $(FULLFLAG) $^ -o $(EXE) $(CFLAGS) -lfftw3 -lstdc++; mv $(EXE) bin/

clean:
	rm -rf *~ ${EXE} *.o *.gch

