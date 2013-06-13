input: clean
	gmsh -2 -algo del2d src/square.geo
	cp src/square.msh .

clean:
	rm -f *.vtu *.msh *.s *.d.1 *.stat \
	matrixdump matrixdump.info
	rm -f src/*.msh
