input: clean
	gmsh -2 ./src/KatoPhillips-2D-periodic.geo
	cp src/KatoPhillips-2D-periodic.msh .

clean: 
	rm -f src/*.msh *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump*
