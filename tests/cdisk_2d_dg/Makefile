input: clean
	gmsh -2 -o cdisk.msh src/cdisk.geo 
	../../bin/gmsh2triangle --2d cdisk.msh

clean:
	rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh
