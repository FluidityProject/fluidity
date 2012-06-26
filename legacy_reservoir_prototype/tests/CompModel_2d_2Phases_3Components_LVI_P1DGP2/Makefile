SHELL = sh

SIM = CompModel_2d_2Phases_3Components_P1DGP2_LVI

input: clean

	gmsh -2 src/2d_mesh.geo
	../../../scripts/gmsh2triangle --2d src/2d_mesh.msh

clean:
	rm -f *.d.* fluidity.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
