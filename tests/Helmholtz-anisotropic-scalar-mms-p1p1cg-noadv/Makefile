input:	clean
	gmsh -2 -algo frontal -o MMS_A.msh src/MMS_A_structured.geo
	gmsh -2 -algo frontal -o MMS_B.msh src/MMS_B_structured.geo
	gmsh -2 -algo frontal -o MMS_C.msh src/MMS_C_structured.geo
	gmsh -2 -algo frontal -o MMS_D.msh src/MMS_D_structured.geo
	gmsh -2 -algo frontal -o MMS_E.msh src/MMS_E_structured.geo
	../../bin/gmsh2triangle --2d MMS_A.msh
	../../bin/gmsh2triangle --2d MMS_B.msh
	../../bin/gmsh2triangle --2d MMS_C.msh
	../../bin/gmsh2triangle --2d MMS_D.msh
	../../bin/gmsh2triangle --2d MMS_E.msh

clean:
	rm -rf *.vtu *.stat *.node *.ele *.edge *checkpoint*
