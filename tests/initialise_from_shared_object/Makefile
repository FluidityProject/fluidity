default: input

input: clean
	gmsh -2 src/box_2d.geo

	gfortran -fPIC -shared -g -o example.so src/example.F90 

clean:
	rm -f *.vtu *.stat src/*.msh *.so
