preprocess:
	@echo **********Creating 1D mesh
	../../bin/interval --dx=100.0 -- 0.0 10000.0 line
	@echo **********Creating pseudo-1D mesh
	gmsh -2 -o enforce_hydrostatic_balance_pseudo1d.msh src/enforce_hydrostatic_balance_pseudo1d.geo
	../../bin/gmsh2triangle --2d enforce_hydrostatic_balance_pseudo1d.msh

run:
	@echo **********Running 1D simulation
	../../bin/fluidity enforce_hydrostatic_balance_1d.flml
	@echo **********Running pseudo-1D simulation
	../../bin/fluidity enforce_hydrostatic_balance_pseudo1d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.bound
	rm -f matrixdump* *.log* *.err*

