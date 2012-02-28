input:	clean
	gmsh -3 src/MMS_A.geo -o MMS_A.msh
	gmsh -3 src/MMS_B.geo -o MMS_B.msh
	gmsh -3 src/MMS_C.geo -o MMS_C.msh
	gmsh -3 src/MMS_D.geo -o MMS_D.msh

clean:
	rm -f *.vtu *.stat *.log *checkpoint* *.msh
