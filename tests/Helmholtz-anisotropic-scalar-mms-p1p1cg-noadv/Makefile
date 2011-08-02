PYTHONPATH=PYTHONPATH:/data/JBull/fluidity/turbulence/python/
input:	clean
	gmsh -2 -algo frontal -o MMS_A.msh src/MMS_A.geo
	gmsh -2 -algo frontal -o MMS_B.msh src/MMS_B.geo
	gmsh -2 -algo frontal -o MMS_C.msh src/MMS_C.geo
	gmsh -2 -algo frontal -o MMS_D.msh src/MMS_D.geo
	gmsh -2 -algo frontal -o MMS_E.msh src/MMS_E.geo
	../../scripts/gmsh2triangle --2d MMS_A.msh
	../../scripts/gmsh2triangle --2d MMS_B.msh
	../../scripts/gmsh2triangle --2d MMS_C.msh
	../../scripts/gmsh2triangle --2d MMS_D.msh
	../../scripts/gmsh2triangle --2d MMS_E.msh

clean:
	rm -rf *.vtu *.stat *.node *.ele *.edge *checkpoint*
