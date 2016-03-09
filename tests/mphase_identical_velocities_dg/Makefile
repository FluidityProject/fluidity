preprocess:
	@echo **********Converting mesh from gmsh to triangle format
	cp src/mphase_identical_velocities_dg.msh .

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_identical_velocities_dg.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* mphase_*.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

