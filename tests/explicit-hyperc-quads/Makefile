input: clean
	gmsh -2 src/channel.geo
	../../bin/gmsh2triangle --2d src/channel.msh

clean:
	rm -rf *.d.* *.s* *.vtu *.edge channel.node channel.ele channel.face *checkpoint* src/channel.msh
