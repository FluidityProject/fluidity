input:	clean
	gmsh2triangle --2d src/MMS_A.msh
	gmsh2triangle --2d src/MMS_B.msh
	gmsh2triangle --2d src/MMS_C.msh
	gmsh2triangle --2d src/MMS_D.msh
	gmsh2triangle --2d src/MMS_E.msh

clean:
	rm -f *.vtu *.stat *.log* *.err* *.node *.ele *.edge *checkpoint* *.convergence matrixdump*
