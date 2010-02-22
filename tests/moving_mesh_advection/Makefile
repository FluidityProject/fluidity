input: clean
	gmsh -2 src/channel.geo
	../../scripts/gmsh2triangle --2d src/channel.msh

clean:
	rm -rf *.stat *.vtu src/channel.msh channel.node channel.ele channel.edge *.log-0 *.err-0
