input: clean
	gmsh -format msh2 -2 src/MMS_B.geo -o MMS_B.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.msh *.ele *.node *.edge
	rm -f *.stat *_periodised.flml \
	matrixdump matrixdump.info
