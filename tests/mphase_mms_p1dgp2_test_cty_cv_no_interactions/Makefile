preprocess:
	@echo **********Creating 2D meshes
	../../bin/gmsh2triangle --2d src/MMS_A.msh
	../../bin/gmsh2triangle --2d src/MMS_B.msh
	../../bin/gmsh2triangle --2d src/MMS_C.msh

run:
	@echo **********Running simulations
	../../bin/fluidity -v3 -l MMS_A.flml
	../../bin/fluidity -v3 -l MMS_B.flml
	../../bin/fluidity -v3 -l MMS_C.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly
	rm -f matrixdump* *checkpoint* *.log* *.err*
