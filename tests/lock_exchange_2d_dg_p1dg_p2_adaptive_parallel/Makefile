input: clean
	gmsh -2 ./src/channel_2d.geo
	../../scripts/gmsh2triangle --2d ./src/channel_2d.msh
	fldecomp -n 2 -i triangle -f channel_2d

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.face *.poly *vtu *.s  *.stat *.log-? *.err-? matrixdump* *checkpoint* *.halo
