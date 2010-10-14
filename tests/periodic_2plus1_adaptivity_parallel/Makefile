input: clean
	gmsh -2 src/MMS_B.geo -o MMS_B.msh
	../../scripts/gmsh2triangle --2d MMS_B.msh
	../../bin/fldecomp -n 3 MMS_B

clean:
	rm -f *.vtu *.log-? *.err-? *.msh *.halo *.ele *.node *.edge
	rm -f *.stat *_periodised.flml
