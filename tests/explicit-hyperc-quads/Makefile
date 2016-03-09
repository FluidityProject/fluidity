input: clean
	gmsh -2 src/channel.geo
	cp src/channel.msh .

clean:
	rm -rf *.d.* *.s* *.vtu *.edge channel.node channel.ele channel.face *checkpoint* src/channel.msh
