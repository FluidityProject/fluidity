PROJECT = channel

input: clean
	gmsh -3 -bin src/channel.geo
	mv src/channel.msh .

clean:
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.msh
