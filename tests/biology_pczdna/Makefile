input: clean
	gmsh -3 -o column-3D.msh src/column.geo
	../../scripts/gmsh2triangle column-3D.msh

clean:
	rm -f  *.ele *.edge *.face *.msh *.node *.poly *.vtu *.s *.d.1 *.stat *.detectors
