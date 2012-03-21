CPPFLAGS = -I./include/ -I/usr/include/
CFLAGS = -Wall -g -std=gnu99
FFLAGS = -Wall -g
AR = ar rcv
RANLIB = ranlib
LIBS = -L./ -L./lib/libexodus.a -L/usr/lib/ -lnetcdf



FC = gfortran
CC = gcc

default: modules libb

libb:
	$(AR) libexodusiifortran.a exodusiicwrapper.o exodusii_fortran.mod
	$(RANLIB) libexodusiifortran.a



modules: exodusii_fortran.o exodusiicwrapper.o
exodusii_fortran.o: exodusii_fortran.F90
	$(FC) -c exodusii_fortran.F90
exodusiicwrapper.o: exodusiicwrapper.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c exodusiicwrapper.c

clean:
	rm -f *.mod *.o *.a

distclean: clean
	rm -f make.log include/exodusII.h include/exodusII_ext.h include/exodusII_int.h
