input: clean
	gmsh -3 src/mmat-flowthrough.geo
	../../scripts/gmsh2triangle src/mmat-flowthrough.msh

clean:
	rm -rf *.d.* *.s* *.vtu mmat-flowthrough.node mmat-flowthrough.ele mmat-flowthrough.face 
