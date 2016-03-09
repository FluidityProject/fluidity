input: clean
	gmsh -2 src/electrochemical.geo
	cp src/electrochemical.msh .

clean:
	rm -rf *.stat *.vtu *.convergence *.node *.edge *.ele *checkpoint* *.log-0 *.err-0 ./src/*.msh matrix* prolon*
