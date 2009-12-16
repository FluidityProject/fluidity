input: clean
	gmsh -3 src/elastic-cube.geo
	../../scripts/gmsh2triangle src/elastic-cube.msh

clean:
	rm -rf *.d.* *.s* *.vtu elastic-cube.fld elastic-cube0 elastic-cube elastic-cube.node elastic-cube.ele elastic-cube.face
