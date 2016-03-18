input: clean
	gmsh -2 src/annulus_2d.geo
	cp src/annulus_2d.msh .

clean:
	rm -f *.vtu *.stat *.msh matrixdump matrixdump.info *.log-0 *.err-0\
	  src/*.msh
