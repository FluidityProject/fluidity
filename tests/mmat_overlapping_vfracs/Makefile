input: clean
	gmsh -2 src/column_2d.geo
	cp src/column_2d.msh .

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint* src/column_2d.msh fluidity.[le]* *.detectors.h5part \
	matrixdump matrixdump.info
