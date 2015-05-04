preprocess:
	@echo **********Creating 2D mesh
	cp src/bc_from_field.msh .

run:
	@echo **********Running simulation
	../../bin/fluidity -v2 -l scalar_field_bc_from_field.flml
	../../bin/fluidity -v2 -l vector_field_bc_from_field.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *.log* *.err*

