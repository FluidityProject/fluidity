input: clean
	gmsh -3 src/mmat-gravity-col.geo
	cp src/mmat-gravity-col.msh .

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence mmat-gravity-col.node mmat-gravity-col.face mmat-gravity-col.ele src/set_phi *checkpoint* src/mmat-gravity-col.msh fluidity.[le]* \
	matrixdump matrixdump.info
