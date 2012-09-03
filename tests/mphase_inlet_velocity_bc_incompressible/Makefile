preprocess:
	@echo **********Creating 1D mesh
	../../bin/interval --dx=0.001 -- 0.0 0.1 line

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l mphase_inlet_velocity_bc_incompressible.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.bound
	rm -f matrixdump* *.log* *.err*

