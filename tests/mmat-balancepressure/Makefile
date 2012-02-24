input: clean
	gmsh -2 src/2d_square.geo
	../../bin/gmsh2triangle --2d src/2d_square.msh

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence *.edge *.node *.face *.ele *checkpoint* matrixdump* src/2d_square.msh *.log 
