input:	clean
	gmsh -format msh2 -2 src/MMS_A.geo -o MMS_A.msh
	gmsh -format msh2 -2 src/MMS_B.geo -o MMS_B.msh
	gmsh -format msh2 -2 src/MMS_C.geo -o MMS_C.msh
	gmsh -format msh2 -2 src/MMS_D.geo -o MMS_D.msh

clean:
	rm -f *.vtu *.stat *.log *checkpoint* *.msh
	rm -f *.err-? *.log-? \
	matrixdump matrixdump.info
