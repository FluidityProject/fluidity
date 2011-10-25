input: clean
	gmsh -2 ./src/channel_2d.geo
	../../bin/gmsh2triangle --2d ./src/channel_2d.msh
	../../bin/fldecomp -n 2 -f channel_2d

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.face *.poly *vtu *.s  *.stat *.log-* *.err-* matrixdump*
