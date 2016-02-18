input: clean
	gmsh -2 src/square.geo -o square.msh

clean:
	rm -rf *.msh *.stat *.vtu *.log-0 *.err-0 *checkpoint *.convergence *.log *.node *.edge *.ele *.halo *.pvtu viscous_fs_zhong_B_* viscous_fs_zhong_A_* fluidity.* *.pyc ib_* flredecomp* parallel*
	rm -rf *checkpoint* matrixdump*
