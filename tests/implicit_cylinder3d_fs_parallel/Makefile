input: clean
	gmsh -2 box2d-coarse.msh src/box2d-coarse.geo
	gmsh -3 cylinder-coarse.msh src/cylinder-coarse.geo
	../../scripts/gmsh2triangle --2d box2d-coarse.msh 
	../../scripts/gmsh2triangle cylinder-coarse.msh 
	../../bin/fldecomp -n 2 box2d-coarse

clean:
	rm -rf *.d* *.s* *.pvtu *.log *.node *.ele *.edge *.face *.vtu matrixdump* *.pyc *.msh fluidity.* *.halo *.dat *_0 *_1 *_2
