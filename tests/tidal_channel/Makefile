PROJECT = channel

input: clean
	gmsh -3 src/channel.geo
	mv src/channel.msh .
	../../scripts/gmsh2triangle channel.msh

clean:
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node *.msh
