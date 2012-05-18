SHELL = sh

SIM = BL_test_10elements_2d

input: clean

	gmsh -2 src/2d_mesh.geo
	../../../scripts/gmsh2triangle --2d src/2d_mesh.msh

clean:
	rm -f *.d.* fluidity.* ocvn
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat matrixdump matrixdump.info
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
