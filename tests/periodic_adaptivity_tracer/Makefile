input:
	gmsh -2 src/mesh.geo -o mesh.msh
	../../scripts/gmsh2triangle --2d mesh.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.stat *.msh
