input: clean
	gmsh -3 src/mmat-shocktube.geo
	cp src/mmat-shocktube.msh .

clean:
	rm -rf *.d.* *.stat *.vtu mmat-shocktube.node mmat-shocktube.ele mmat-shocktube.face *checkpoint* *.convergence fluidity.[le]* src/mmat-shocktube.msh \
	matrixdump matrixdump.info
