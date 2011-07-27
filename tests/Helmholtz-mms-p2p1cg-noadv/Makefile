PYTHONPATH=PYTHONPATH:/data/JBull/fluidity/turbulence/python/
input:	clean
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo
	gmsh -2 src/MMS_D.geo
	../../scripts/gmsh2triangle --2d src/MMS_A.msh
	../../scripts/gmsh2triangle --2d src/MMS_B.msh
	../../scripts/gmsh2triangle --2d src/MMS_C.msh
	../../scripts/gmsh2triangle --2d src/MMS_D.msh

clean:
	rm -rf *.vtu *.stat *.node *.ele *.edge *checkpoint*
