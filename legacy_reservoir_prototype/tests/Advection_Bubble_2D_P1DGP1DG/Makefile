SHELL = sh

SIM = Advection_Bubble_P1DGP1

input: clean

	gmsh -2 src/2d_mesh.geo
	../../../scripts/gmsh2triangle --2d src/2d_mesh.msh

clean:
	rm -f *.d.* fluidity.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.detectors *.vtu
	rm -f mirror_new.dat matrixdump*
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
