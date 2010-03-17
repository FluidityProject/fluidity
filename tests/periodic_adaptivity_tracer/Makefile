input:
	gmsh -2 src/mesh.geo
	../../scripts/gmsh2triangle --2d src/mesh.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.stat
