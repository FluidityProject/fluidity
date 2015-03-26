input: clean
	@echo "\n********** Building surface mesh with Gmsh :"
	gmsh -2 -o sphericalSegment.msh src/sphericalSegment.geo
	@echo "\n********** Decomposing the mesh into 2 parts for parallel run:"
	mpirun -n 2 ../../bin/flredecomp -i 1 -o 2 sphericalSegment sphericalSegment_flredecomp

clean:
	rm -f  *.vtu *.stat *.msh *_flredecomp* matrixdump*
