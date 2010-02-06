PROJECT = channel

input: clean
	gmsh -3 -o channel.msh src/channel.geo
	../../scripts/gmsh2triangle channel.msh
	../../bin/periodise channel.flml

clean:
	rm -f $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node 
	rm -f *.msh
	rm -f *_periodised.flml
