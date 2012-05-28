input: clean
	gmsh -3 -o column.msh src/column.geo
	../../bin/gmsh2triangle column.msh
	../../bin/fldecomp -n 3 column

clean: 
	rm -rf *.halo *.vtu *.pvtu *.log* *.err* *.stat *.detectors *.dat *.ele *.node *.bound
	rm -rf vew_toymodel_*
