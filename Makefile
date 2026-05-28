# make file for Kenni Dinese Petersen's code

# if on glmac remember to run >> source /opt/intel/bin/compilervars.sh intel64
# if on malbec remember to run >> module load intel/compiler
# if on local remember to run >> module load GCCcore/12.3.0

## glmac options
#CC = icc
#CFLAGS = -O2 -xHOST -openmp

## malbec options
#CC = icc
#CFLAGS = -O2 -Wall -xAVX -openmp
#LIB = 

## local options
CC = gcc
CFLAGS = -O2 -Wall -fopenmp
LIB = -lm

## debug
#CC = gcc
#CFLAGS = -g -fopenmp
#LIB = -lm

melting:
	$(CC) $(CFLAGS) $(LIB) -o melting.$(CC) main.c \
	momentummarkerfuncs_par.c multigridfuncs2.c tempfuncs_par.c surfacefuncs.c composition.c
veloforce:
	$(CC) $(CFLAGS) $(LIB) -o veloforce.$(CC) main.c \
	momentummarkerfuncs_par.c multigridfuncs2.c tempfuncs_par.c surfacefuncs.c composition.c
clean:
	rm melting.$(CC) veloforce.$(CC)