input:	clean
	../../bin/gmsh2triangle -2 MMS_A_unstructured.msh
	../../bin/gmsh2triangle -2 MMS_B_unstructured.msh
	../../bin/gmsh2triangle -2 MMS_C_unstructured.msh
	../../bin/gmsh2triangle -2 MMS_D_unstructured.msh

clean:
	rm -f *.vtu *.ele *.node *.edge *.stat *.log fluidity.err-0 fluidity.log-0

flmls:
	@$(foreach M,B C D, \
		cp MMS_A_unstructured.flml MMS_$(M)_unstructured.flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_unstructured.flml > tempfile; \
		mv tempfile MMS_$(M)_unstructured.flml; \
	)
