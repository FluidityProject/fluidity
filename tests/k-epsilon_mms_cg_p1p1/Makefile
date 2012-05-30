input:	clean
	gmsh -2 -bin src/MMS_A.geo
	gmsh -2 -bin src/MMS_B.geo
	gmsh -2 -bin src/MMS_C.geo
	gmsh -2 -bin src/MMS_D.geo

clean:
	rm -rf *.log* *.err* *.vtu *.stat *.node *.ele *.edge *checkpoint* src/*.msh
