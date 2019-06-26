input: clean
	gmsh -format msh2 -2 -o cdisk.msh src/cdisk.geo 

clean:
	rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh
	rm -f matrixdump* *.log-? *.err-?
