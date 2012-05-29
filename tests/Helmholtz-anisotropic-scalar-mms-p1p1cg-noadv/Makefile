input:	clean
	gmsh -2 -algo frontal -o MMS_A_unstructured.msh src/MMS_A_unstructured.geo
	gmsh -2 -algo frontal -o MMS_B_unstructured.msh src/MMS_B_unstructured.geo
	gmsh -2 -algo frontal -o MMS_C_unstructured.msh src/MMS_C_unstructured.geo
	gmsh -2 -algo frontal -o MMS_D_unstructured.msh src/MMS_D_unstructured.geo
	gmsh -2 -algo frontal -o MMS_E_unstructured.msh src/MMS_E_unstructured.geo

clean:
	#rm -rf *.vtu *.stat *.node *.ele *.edge

flmls:
	@$(foreach M,B C D E, \
		cp MMS_A_unstructured.flml MMS_$(M)_unstructured.flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_unstructured.flml > tempfile; \
		mv tempfile MMS_$(M)_unstructured.flml; \
	)
