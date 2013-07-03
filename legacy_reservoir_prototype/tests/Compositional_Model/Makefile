SHELL = sh

SIM = Buckley_Leverett_2D

input: clean

	gmsh -2 src/2d_mesh_10layers.geo
	gmsh -2 src/2d_mesh_50layers.geo
	gmsh -2 src/2d_mesh_100layers.geo
	gmsh -2 src/2d_mesh_250layers.geo
	gmsh -2 src/2d_mesh_500layers.geo

	../../../scripts/gmsh2triangle --2d src/2d_mesh_10layers.msh
	../../../scripts/gmsh2triangle --2d src/2d_mesh_50layers.msh
	../../../scripts/gmsh2triangle --2d src/2d_mesh_100layers.msh
	../../../scripts/gmsh2triangle --2d src/2d_mesh_250layers.msh
	../../../scripts/gmsh2triangle --2d src/2d_mesh_500layers.msh

clean:
	rm -f *.d.* fluidity.*
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat matrixdump*
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
	rm -f *_50layers.mpml *_100layers.mpml *_250layers.mpml *_500layers.mpml
