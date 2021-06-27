input: clean
	gmsh -3 ./src/channel_3d.geo
	cp src/channel_3d.msh .

clean: 
	rm -f  *.msh *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* src/*.msh
