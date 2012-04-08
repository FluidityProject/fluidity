SHELL = sh

SIM = BL_3d

input: clean

	gmsh -optimize -3 src/3d_mesh.geo
	../../../scripts/gmsh2triangle src/3d_mesh.msh

clean:
	rm -f *.d.* fluidity.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
