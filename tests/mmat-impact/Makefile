input: clean
	gmsh -2 src/mmat-impact.geo
	../../scripts/gmsh2triangle --2d src/mmat-impact.msh

clean:
	rm -rf *.d.* *.s* 2material_impact_*.vtu mmat-impact.node mmat-impact.ele mmat-impact.face *.convergence *checkpoint* output matrixdump*
