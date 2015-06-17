preprocess:
	@echo **********Creating 2D meshes
	cp src/MMS_A.msh .
	cp src/MMS_B.msh .
	cp src/MMS_C.msh .
	cp src/MMS_D.msh .

run:
	@echo **********Running simulations
	../../bin/fluidity -v3 -l MMS_A.flml
	../../bin/fluidity -v3 -l MMS_B.flml
	../../bin/fluidity -v3 -l MMS_C.flml
	../../bin/fluidity -v3 -l MMS_D.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *checkpoint* *.log* *.err*
