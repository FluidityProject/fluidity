input: clean
	gmsh -2 src/2d_square.geo
	cp src/2d_square.msh .

clean:
	rm -rf *.stat *.vtu *.node *.ele *.edge *checkpoint* *convergence* *.log-0 *.err-0 src/2d_square.msh \
	matrixdump matrixdump.info
