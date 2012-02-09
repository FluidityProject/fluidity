input:	clean
	gmsh -3 src/box.geo
	../../bin/gmsh2triangle src/box.msh

clean:
	rm -f *.vtu *.stat *.s *.log* *.err* src/*.msh *.node *.ele *.edge *.face *checkpoint* matrixdump*
