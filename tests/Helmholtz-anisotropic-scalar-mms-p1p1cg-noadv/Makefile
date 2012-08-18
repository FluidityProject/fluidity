MESH=structured
input:	clean
	gmsh -2 -optimize -algo frontal -o MMS_A_$(MESH).msh src/MMS_A_$(MESH).geo
	gmsh -2 -optimize -algo frontal -o MMS_B_$(MESH).msh src/MMS_B_$(MESH).geo
	gmsh -2 -optimize -algo frontal -o MMS_C_$(MESH).msh src/MMS_C_$(MESH).geo
	gmsh -2 -optimize -algo frontal -o MMS_D_$(MESH).msh src/MMS_D_$(MESH).geo
	gmsh -2 -optimize -algo frontal -o MMS_E_$(MESH).msh src/MMS_E_$(MESH).geo

clean:
	rm *.msh

flmls:
	@$(foreach M,B C D E, \
		cp MMS_A_$(MESH).flml MMS_$(M)_$(MESH).flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_$(MESH).flml > tempfile; \
		mv tempfile MMS_$(M)_$(MESH).flml; \
	)
