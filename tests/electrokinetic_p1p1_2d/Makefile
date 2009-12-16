input: clean
	gmsh -2 src/electrokinetic.geo
	../../scripts/gmsh2triangle --2d src/electrokinetic.msh

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint* *.log-0 *.err-0 ./src/*.msh matrix* prolon*
