preprocess:
	@echo **********Creating 2D meshes
	cp src/MMS_A.msh .
	cp src/MMS_B.msh .
	cp src/MMS_C.msh .

run:
	@echo **********Running simulations
	mpiexec fluidity MMS_A.flml
	mpiexec fluidity MMS_B.flml
	mpiexec fluidity MMS_C.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.poly
	rm -f matrixdump* *checkpoint* *.log* *.err* *.convergence
