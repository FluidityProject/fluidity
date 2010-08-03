input: clean
	gmsh -3 src/channel.geo
	../../scripts/gmsh2triangle src/channel.msh

clean:
	rm -rf *.stat *.vtu channel.node channel.ele channel.face *.convergence *checkpoint* fluidity.log-0 fluidity.err-0
