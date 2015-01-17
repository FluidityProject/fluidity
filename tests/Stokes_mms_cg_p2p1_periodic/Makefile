input:	clean
	gmsh -2 src/MMS_A.geo -o MMS_A.msh
	gmsh -2 src/MMS_B.geo -o MMS_B.msh
	gmsh -2 src/MMS_C.geo -o MMS_C.msh
	gmsh -2 src/MMS_D.geo -o MMS_D.msh

clean:
	rm -f *.vtu *.stat *.log *checkpoint* *.msh
	rm -f *.err-? *.log-? \
	matrixdump matrixdump.info
