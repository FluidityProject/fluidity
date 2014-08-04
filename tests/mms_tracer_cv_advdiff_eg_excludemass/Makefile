input:	clean
	gmsh2triangle --2d src/MMS_A.msh
	gmsh2triangle --2d src/MMS_B.msh
	gmsh2triangle --2d src/MMS_C.msh
	gmsh2triangle --2d src/MMS_D.msh
	gmsh2triangle --2d src/MMS_E.msh
	#gmsh2triangle --2d src/MMS_F.msh

clean:
	rm -f *.vtu *.stat *.log *.node *.ele *.edge *checkpoint* *.convergence \
	matrixdump matrixdump.info
