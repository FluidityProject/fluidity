input: default
default: clean
	gmsh -2 square.geo
	../../scripts/gmsh2triangle --2d square.msh
	fldecomp -n 2 square
clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat