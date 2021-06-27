input: clean
	gmsh -2 src/drainage_a.geo
	gmsh -2 src/drainage_b.geo
	gmsh -2 src/drainage_c.geo
	gmsh -2 src/drainage_d.geo
	cp src/drainage_a.msh .
	cp src/drainage_b.msh .
	cp src/drainage_c.msh .
	cp src/drainage_d.msh .

clean:
	rm -rf *.stat *.vtu *.msh *.log-0 *.err-0
	rm -rf matrixdump*
