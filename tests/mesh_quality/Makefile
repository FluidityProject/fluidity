input: clean
	gmsh -format msh2 -2 src/square.geo .

clean:
	rm -f *.ele *.edge *.node src/*.msh *.poly *.vtu *.s *.d.1 *.stat	fluidity.[le]* \
	matrixdump matrixdump.info
