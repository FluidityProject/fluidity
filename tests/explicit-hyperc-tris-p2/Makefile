input: clean
	gmsh -2 src/channel.geo
	cp src/channel.msh .

clean:
	rm -rf *.d.* *.s* *.vtu channel.node channel.ele channel.edge *checkpoint* *.convergence src/*.msh *.log-0 *.err-0 *.detectors \
	matrixdump matrixdump.info
