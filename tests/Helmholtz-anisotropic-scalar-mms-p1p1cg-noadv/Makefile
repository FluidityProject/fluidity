input:	clean
	gmsh -2 -algo frontal -o MMS_A.msh src/MMS_A.geo
	gmsh -2 -algo frontal -o MMS_B.msh src/MMS_B.geo
	gmsh -2 -algo frontal -o MMS_C.msh src/MMS_C.geo
	gmsh -2 -algo frontal -o MMS_D.msh src/MMS_D.geo
	../../bin/gmsh2triangle --2d MMS_A.msh
	../../bin/gmsh2triangle --2d MMS_B.msh
	../../bin/gmsh2triangle --2d MMS_C.msh
	../../bin/gmsh2triangle --2d MMS_D.msh

clean:
	rm -rf *.vtu *.stat *.node *.ele *.edge *checkpoint*
