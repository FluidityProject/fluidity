preprocess:
	@echo **********Creating 2D meshes
	../../bin/gmsh2triangle --2d src/MMS_A.msh
	../../bin/gmsh2triangle --2d src/MMS_B.msh
	../../bin/gmsh2triangle --2d src/MMS_C.msh
	../../bin/gmsh2triangle --2d src/MMS_D.msh

run:
	@echo **********Running simulations
	mpiexec ../../bin/fluidity MMS_A.flml
	mpiexec ../../bin/fluidity MMS_B.flml
	mpiexec ../../bin/fluidity MMS_C.flml
	mpiexec ../../bin/fluidity MMS_D.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *checkpoint* *.log* *.err* *.convergence
