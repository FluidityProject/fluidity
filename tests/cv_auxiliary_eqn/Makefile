input: clean
	#gmsh -2 src/MMS_A.geo
	#gmsh -2 src/MMS_B.geo
	#gmsh -2 src/MMS_C.geo
	#gmsh -2 src/MMS_D.geo
	#gmsh -2 src/MMS_E.geo
	#gmsh -2 src/MMS_F.geo
	../../scripts/gmsh2triangle --2d src/MMS_A.msh
	../../scripts/gmsh2triangle --2d src/MMS_B.msh
	../../scripts/gmsh2triangle --2d src/MMS_C.msh
	../../scripts/gmsh2triangle --2d src/MMS_D.msh
	../../scripts/gmsh2triangle --2d src/MMS_E.msh
	../../scripts/gmsh2triangle --2d src/MMS_F.msh

clean:
	rm -rf *.stat *.vtu *.node *.ele *.face *.convergence *checkpoint*
