PROJECT = rt0_projection_quads

input: clean
	gmsh -2 box.geo
	../../tools/gmsh2triangle.py box.msh

clean:	
	rm -f *.vtu
	rm -f *.ele *.face *.node *.edge *.stat
	rm -f *.msh
	rm -f log