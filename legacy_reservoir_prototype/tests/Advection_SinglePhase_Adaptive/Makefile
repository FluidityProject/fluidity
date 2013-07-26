SHELL = sh

#SIM = Advection_SinglePhase_Adaptive

input: clean

	gmsh -2 src/2d_mesh_unstructured.geo
	gmsh -2 src/2d_mesh_adaptive.geo
	gmsh -2 src/2d_mesh_square.geo
	../../../scripts/gmsh2triangle --2d src/2d_mesh_adaptive.msh
	../../../scripts/gmsh2triangle --2d src/2d_mesh_unstructured.msh
	../../../scripts/gmsh2triangle --2d src/2d_mesh_square.msh

clean:
	rm -f *.d.* fluidity.* ocvn
	rm -f flog.dat time_count tempfile CheckSparsityMatrix.dat filedim fvalues mirror_int_data.dat fort.*
	rm -f *.node *.bound *.ele *.face *.stat *.vtu
	rm -f mirror_new.dat matrixdump matrixdump.info
	rm -f *.ele *.bound *.face *.node *.edge
	rm -f src/*.msh
