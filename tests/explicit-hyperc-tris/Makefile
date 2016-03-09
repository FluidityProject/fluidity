input: clean
	gmsh -2 src/channel.geo
	cp src/channel.msh .

clean:
	rm -rf *.d.* *.s* *.vtu channel.node channel.ele channel.edge *.log-0 *.err-0 *checkpoint* src/channel.msh Water__Tracer.convergence \
	matrixdump matrixdump.info
