preprocess:
	@echo **********Creating 1D mesh
	../../bin/interval --dx=0.00025 -- -0.3 0.3 line

run:
	@echo **********Running simulation
	../../bin/fluidity mphase_dusty_gas_shock_tube.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.bound
	rm -f matrixdump* *.log* *.err*

