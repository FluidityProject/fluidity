input: clean
	gmsh -3 box3d-coarse.msh src/box3d-coarse.geo
	gmsh -3 cylinder-coarse.msh src/cylinder-coarse.geo
	../../scripts/gmsh2triangle box3d-coarse.msh 
	../../scripts/gmsh2triangle cylinder-coarse.msh 
	../../bin/fldecomp -n 4 box3d-coarse 

clean:
	rm -rf *.d* *.s* *.pvtu *.log *.node *.ele *.edge *.face *.vtu matrixdump* *.pyc *.msh fluidity.* *.halo *.dat *_0 *_1 *_2 *_3
