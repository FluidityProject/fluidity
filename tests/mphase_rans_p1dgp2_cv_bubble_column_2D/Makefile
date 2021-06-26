preprocess:
	@echo **********Creating 2D mesh
	gmsh -2 -o bubcol2d_pfleger.msh src/bubcol2d_pfleger.geo

run:
	@echo **********Running flredecomp to decompose the domain
	mpiexec -n 4 flredecomp -i 1 -o 4 mphase_rans_p1dgp2_cv_bubble_column_2D mphase_rans_p1dgp2_cv_bubble_column_2D_flredecomp
	@echo **********Running simulation
	mpiexec -n 4 fluidity -v2 -l mphase_rans_p1dgp2_cv_bubble_column_2D_flredecomp.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.vtu
	rm -f *.msh
	rm -f matrixdump* *.log* *.err*
	rm -rf bubcol_2d_rans_*
	rm -f *flredecomp*

