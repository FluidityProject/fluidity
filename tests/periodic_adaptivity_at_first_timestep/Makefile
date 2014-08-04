input:
	gmsh -2 src/MMS_B.geo -o MMS_B.msh
	gmsh2triangle --2d MMS_B.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.msh *.edge *.ele *.node *.stat *_periodised.flml
	rm -f matrixdump*
