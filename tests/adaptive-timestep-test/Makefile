input: clean
	gmsh -3 src/channel.geo
	cp src/channel.msh .

clean:
	rm -rf *.d.* *.s* *.vtu channel.node channel.ele channel.face *.convergence *checkpoint* src/channel.msh \
	matrixdump matrixdump.info
