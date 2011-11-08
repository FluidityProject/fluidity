PROJECT = channel

input: clean
	gmsh -2 src/channel.geo
	mv src/channel.msh .
	../../bin/gmsh2triangle --2d channel.msh

clean:
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.edge *.node *.msh
