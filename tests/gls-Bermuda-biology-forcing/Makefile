input: clean
	gmsh -3 -o gls-Bermuda.msh src/gls-Bermuda.geo
	../../scripts/gmsh2triangle gls-Bermuda.msh

clean: 
	rm -f  *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *log *err
