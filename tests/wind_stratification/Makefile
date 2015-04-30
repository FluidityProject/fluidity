input: clean
	gmsh -2 src/square.geo -o square.msh
	../../bin/transform_mesh '(x,0.001*y)' square

clean:
	rm -f  *.ele *.edge *.node *.msh *.vtu *.stat matrixdump* *.node.bak fluidity.*
