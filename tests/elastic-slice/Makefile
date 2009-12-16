input: clean
	gmsh -3 src/elastic-slice.geo
	../../scripts/gmsh2triangle src/elastic-slice.msh

clean:
	rm -rf *.d.* *.s* *.vtu elastic-slice.node elastic-slice.ele elastic-slice.face
