input: clean
	gmsh -2 src/square.geo -o square.msh
	../../scripts/gmsh2triangle --2d square.msh
	../../scripts/transform_mesh '(x,0.001*y)' square

clean:
	rm -f  *.ele *.edge *.node *.msh *.vtu *.stat matrixdump* *.node.bak fluidity.*
