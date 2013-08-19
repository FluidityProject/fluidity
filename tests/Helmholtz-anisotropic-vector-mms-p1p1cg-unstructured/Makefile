input:	clean
	gmsh -2 -bin -algo front2d -o MMS_A_unstructured.msh src/MMS_A_unstructured.geo
	gmsh -2 -bin -algo front2d -o MMS_B_unstructured.msh src/MMS_B_unstructured.geo
	gmsh -2 -bin -algo front2d -o MMS_C_unstructured.msh src/MMS_C_unstructured.geo
	gmsh -2 -bin -algo front2d -o MMS_D_unstructured.msh src/MMS_D_unstructured.geo

clean:
	rm -f *.vtu *.msh *.stat *.log fluidity.log-0 fluidity.err-0

flmls:
	@$(foreach M,B C D, \
		cp MMS_A_unstructured.flml MMS_$(M)_unstructured.flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_unstructured.flml > tempfile; \
		mv tempfile MMS_$(M)_unstructured.flml; \
	)
