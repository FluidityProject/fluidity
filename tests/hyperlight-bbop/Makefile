input: clean
	gmsh -3 -o hyperlight-bbop.msh src/hyperlight.geo
	../../scripts/gmsh2triangle hyperlight-bbop.msh

clean: 
	rm -f  *.msh *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *log *err
