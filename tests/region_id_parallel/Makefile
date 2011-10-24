input: clean
	gmsh -3 src/region_id.geo
	../../bin/gmsh2triangle src/region_id.msh
	../../bin/fldecomp -n 2  region_id

clean:
	rm -rf *.d.* *.s* *.vtu region_id.node region_id.ele region_id.face
