include ../tools.mk

default: input

input: clean
	gmsh -2 src/box_2d.geo

	$(FC) $(FCFLAGS) -fPIC -shared -o example.so src/example.F90 
	$(CC) $(CFLAGS) -fPIC -shared -o c_example.so src/c_example.c

clean:
	rm -f *.vtu *.stat src/*.msh *.so
