input: clean
	gmsh -2 src/region_id.geo
	../../scripts/gmsh2triangle --2d src/region_id.msh

clean:
	rm -rf *.d.* *.s* *.vtu region_id.node region_id.ele region_id.edge
