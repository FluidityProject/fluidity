input: clean
	gmsh -2 src/channel.geo
	../../scripts/gmsh2triangle --2d src/channel.msh

clean:
	rm -rf *.d.* *.s* *.vtu channel.node channel.ele channel.face
