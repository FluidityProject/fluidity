input:	clean
	gmsh -2 -bin -o MMS_A_structured.msh src/MMS_A_structured.geo
	gmsh -2 -bin -o MMS_B_structured.msh src/MMS_B_structured.geo
	gmsh -2 -bin -o MMS_C_structured.msh src/MMS_C_structured.geo
	gmsh -2 -bin -o MMS_D_structured.msh src/MMS_D_structured.geo

clean:
	rm -f *.vtu *.msh *.stat *.log fluidity.err-0 fluidity.log-0

flmls:
	@$(foreach M,B C D, \
		cp MMS_A_structured.flml MMS_$(M)_structured.flml; \
		sed 's/MMS_A/MMS_'$(M)'/' MMS_$(M)_structured.flml > tempfile; \
		mv tempfile MMS_$(M)_structured.flml; \
	)
