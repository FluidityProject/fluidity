input: clean
	gmsh -2 -o square.msh src/square.geo && ../../bin/gmsh2triangle --2d square.msh

clean:
	rm -f *.msh *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat matrixdump* *log* *err*
