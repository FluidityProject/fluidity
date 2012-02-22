PYTHONPATH=PYTHONPATH:/data/JBull/fluidity/turbulence/python/
input:	clean
	gmsh -2  -o MMS_A.msh src/MMS_A.geo
	gmsh -2  -o MMS_B.msh src/MMS_B.geo
	gmsh -2  -o MMS_C.msh src/MMS_C.geo
	gmsh -2  -o MMS_D.msh src/MMS_D.geo
	gmsh -2  -o MMS_E.msh src/MMS_E.geo
	#gmsh -2  -o MMS_B.msh src/MMS_B_test.geo
	#gmsh -2 -algo frontal -o MMS_A.msh src/MMS_A.geo
	../../bin/gmsh2triangle --2d MMS_A.msh
	../../bin/gmsh2triangle --2d MMS_B.msh
	../../bin/gmsh2triangle --2d MMS_C.msh
	../../bin/gmsh2triangle --2d MMS_D.msh
	../../bin/gmsh2triangle --2d MMS_E.msh

clean:
	#rm -rf *.vtu *.stat *.node *.ele *.edge *checkpoint*

flmls:
	@$(foreach M,B C D E, \
		cp MMS_A_Helmiso.flml MMS_$(M)_Helmiso.flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_Helmiso.flml > tempfile; \
		mv tempfile MMS_$(M)_Helmiso.flml; \
	)
