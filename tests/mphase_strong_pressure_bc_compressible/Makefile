preprocess:
	@echo **********Creating 1D mesh
	../../bin/interval --dx=100.0 -- 0.0 10000.0 line

run:
	@echo **********Running P2-P1 simulation
	../../bin/fluidity mphase_strong_pressure_bc_compressible_p2p1.flml
	@echo **********Running P0-P1 simulation
	../../bin/fluidity mphase_strong_pressure_bc_compressible_p0p1.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.bound
	rm -f matrixdump* *.log* *.err*

