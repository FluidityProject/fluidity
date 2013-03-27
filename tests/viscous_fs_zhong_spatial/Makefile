input: clean
	gmsh -2 src/squareA.geo
	../../tools/gmsh2triangle.py --2d -i src/squareA.msh
	gmsh -2 src/squareB.geo
	../../tools/gmsh2triangle.py --2d -i src/squareB.msh

clean:
	rm -rf src/*.msh *.stat *.vtu *.log-0 *.err-0 *checkpoint *.convergence *.log *.node *.edge *.ele *.halo *.pvtu viscous_fs_zhong_B_* viscous_fs_zhong_A_* fluidity.* *.pyc
