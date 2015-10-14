#!/bin/csh
#
#
#

CC = icc
COPTIONS = -O3 -xHOST -ipo -openmp -std=c99 -D_GNU_SOURCE
CFLAGS = $(COPTIONS) -I/software/libraries/netcdf/4.1.3-intel/include \
-I/software/libraries/FFTW-3.3/mvapich2-intel

LFLAGS = -L/software/libraries/netcdf/4.1.3-intel/lib -L/software/libraries/FFTW-3.3/mvapich2-intel/lib
LIBS = -lfftw3 -lfftw3_threads -lnetcdf

COMMUN = io_utils.o pml.o propagate.o src.o
OBJECTS =  pve_iso_pml.o $(COMMUN)
OBJECTS2 = pve_vti_pml.o VqP1.o VqP2.o VqS.o $(COMMUN)

bin/pve_iso_pml : $(OBJECTS)
	$(CC) $(CFLAGS) $(LFLAGS) $(LIBS) $(OBJECTS) -o bin/pve_iso_pml

bin/pve_vti_pml : $(OBJECTS2)
	$(CC) $(CFLAGS) $(LFLAGS) $(LIBS) $(OBJECTS2) -o bin/pve_vti_pml

bin/ve_vti_sh_pml : $(COMMUN) ve_vti_sh_pml.o
	$(CC) $(CFLAGS) $(LFLAGS) $(LIBS) $(COMMUN) ve_vti_sh_pml.o -o bin/ve_vti_sh_pml

bin/ve_vti_pml : $(COMMUN) ve_vti_pml.o
	$(CC) $(CFLAGS) $(LFLAGS) $(LIBS) $(COMMUN) ve_vti_pml.o -o bin/ve_vti_pml

io_utils.o : io_utils.c io_utils.h

pml.o : pml.h

propagate.o : propagate.c propagate.h

pve_iso_pml.o : io_utils.h pml.h propagate.h structs.h src.h pve_iso_pml.c

pve_vti_pml.o : io_utils.h pml.h propagate.h structs.h src.h pve_vti_pml.c V.h

src.o : src.c src.h

VqP1.o : V.h VqP1.c

VqP2.o : V.h VqP2.c

VqS.o : V.h VqS.c

ve_vti_pml.o : io_utils.h pml.h propagate.h structs.h src.h  ve_vti_pml.c

ve_vti_sh_pml.o : io_utils.h pml.h propagate.h structs.h src.h  ve_vti_sh_pml.c

clean:
	-rm -f *.o *~ core
