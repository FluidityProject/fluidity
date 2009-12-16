input: clean
	gmsh -2 src/channel.geo
	../../scripts/gmsh2triangle --2d src/channel.msh

clean:
	rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat matrixdump*
