input: clean
	gmsh -3 src/channel.geo
	../../bin/gmsh2triangle src/channel.msh

clean:
	rm -rf *.d.* *.s* *.vtu channel.node channel.ele channel.face *.convergence *checkpoint* src/channel.msh
