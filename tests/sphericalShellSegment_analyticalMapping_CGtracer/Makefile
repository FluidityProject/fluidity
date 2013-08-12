input: clean
	@echo "\n********** Building surface mesh with Gmsh :"
	gmsh -2 -o sphericalSegment.msh src/sphericalSegment.geo
	@echo "\n********** Converting the gmsh mesh (ASCII) to triangle format, mesh is on-the-sphere :"
	../../bin/gmsh2triangle -s sphericalSegment.msh
	@echo "\n********** Decomposing the mesh into 2 parts for parallel run:"
	mpirun -n 2 ../../bin/flredecomp -i 1 -o 2 sphericalSegment sphericalSegment_flredecomp

clean:
	rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh *_flredecomp*
