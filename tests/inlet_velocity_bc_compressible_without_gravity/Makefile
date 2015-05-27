preprocess:
	@echo **********Creating 1D mesh
	../../bin/interval --dx=100.0 -- 0.0 10000.0 line
	@echo **********Creating 2D mesh
	gmsh -2 -o inlet_velocity_bc_compressible_without_gravity_pseudo1d.msh src/inlet_velocity_bc_compressible_without_gravity_pseudo1d.geo

run:
	@echo **********Running simulation
	../../bin/fluidity inlet_velocity_bc_compressible_without_gravity_1d.flml
	../../bin/fluidity inlet_velocity_bc_compressible_without_gravity_pseudo1d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.bound
	rm -f matrixdump* *.log* *.err*

