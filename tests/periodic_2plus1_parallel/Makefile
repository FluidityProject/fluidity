input: clean
	gmsh -2 src/MMS_B.geo -o MMS_B.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.msh *.ele *.node *.edge
	rm -f *.stat *_periodised.flml
	rm -f *vtu
	rm -f *periodised* *flredecomp *2D.flml
	rm -rf periodic_2plus1_parallel_*
