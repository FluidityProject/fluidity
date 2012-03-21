CPPFLAGS = -I./ -I/usr/include/
CFLAGS = -Wall -g -std=gnu99
FFLAGS = -Wall -g
LIBS = -L./ -lexodus -L/usr/lib/ -lnetcdf


FC = gfortran
CC = gcc

default: modules
	$(FC) $(CPPFLAGS) $(FFLAGS) fexowrapper.o cexowrapper.o openexofile.F90 -o openexofile $(LIBS)

modules: fexowrapper.o cexowrapper.o
fexowrapper.o: fexowrapper.F90
	$(FC) -c fexowrapper.F90
cexowrapper.o: cexowrapper.c
	$(CC) $(CFLAGS) -c cexowrapper.c

clean:
	rm -rf *.mod *.o openexofile c_openexofile
