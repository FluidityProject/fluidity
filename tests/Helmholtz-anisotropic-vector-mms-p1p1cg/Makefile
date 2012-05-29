PYTHONPATH=PYTHONPATH:/data/JBull/fluidity/turbulence/python/
input:	clean
	gmsh -2  -o MMS_A_structured.msh src/MMS_A_structured.geo
	gmsh -2  -o MMS_B_structured.msh src/MMS_B_structured.geo
	gmsh -2  -o MMS_C_structured.msh src/MMS_C_structured.geo
	gmsh -2  -o MMS_D_structured.msh src/MMS_D_structured.geo
	gmsh -2  -o MMS_E_structured.msh src/MMS_E_structured.geo
	#gmsh -2 -algo frontal -o MMS_A_structured.msh src/MMS_A_structured.geo

clean:
	#rm -rf *.vtu *.stat *.node *.ele *.edge *checkpoint*

flmls:
	@$(foreach M,B C D E, \
		cp MMS_A_Helmiso_structured.flml MMS_$(M)_Helmiso_structured.flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_Helmiso_structured.flml > tempfile; \
		mv tempfile MMS_$(M)_Helmiso_structured.flml; \
	)
