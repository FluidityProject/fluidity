input: clean
	gmsh -2 src/mmat-impact.geo
	cp src/mmat-impact.msh .

clean:
	rm -rf *.d.* *.s* 2material_impact_*.vtu mmat-impact.node mmat-impact.ele mmat-impact.face *.convergence *checkpoint* output matrixdump* src/mmat-impact.msh fluidity.[le]* mmat-impact.edge
