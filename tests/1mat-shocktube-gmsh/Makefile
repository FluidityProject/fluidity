input: clean
	gmsh -3 -bin src/mmat-shocktube.geo -o mmat-shocktube.msh


clean:
	rm -rf *.d.* *.stat *.vtu mmat-shocktube.node mmat-shocktube.ele mmat-shocktube.face *checkpoint* *.convergence fluidity.[le]* mmat-shocktube.msh \
	matrixdump matrixdump.info
