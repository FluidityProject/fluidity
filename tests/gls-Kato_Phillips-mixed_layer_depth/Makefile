input: clean
	##gmsh -2 ./src/gls-KatoPhillips-2D-periodic.geo
	../../scripts/gmsh2triangle --2d ./src/gls-KatoPhillips-2D-periodic.msh

clean: 
	rm -f  *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump*
