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
	$(AR) libexodusiifortran.a ExodusII_C_Interface.o exodusii_fortran.mod
	$(RANLIB) libexodusiifortran.a


modules: ExodusII_F_Interface.o ExodusII_C_Interface.o
ExodusII_F_Interface.o: ExodusII_F_Interface.F90
	$(FC) -c ExodusII_F_Interface.F90
ExodusII_C_Interface.o: ExodusII_C_Interface.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c ExodusII_C_Interface.c

clean:
	rm -f *.mod *.o *.a

distclean: clean
	rm -f make.log include/exodusII.h include/exodusII_ext.h include/exodusII_int.h
