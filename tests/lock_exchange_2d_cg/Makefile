input: clean
	gmsh -2 ./src/channel_2d.geo
	cp src/channel_2d.msh .

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump*
