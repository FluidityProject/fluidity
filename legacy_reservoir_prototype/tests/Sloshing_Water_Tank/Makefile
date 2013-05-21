SHELL = sh

SIM = Sloshing_Water_Tank

input:  clean

	gmsh -2 src/2d_mesh.geo


clean:
	rm -f *.d.* fluidity.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat matrixdump*
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
	rm -f *_20layers.mpml *_50layers.mpml
