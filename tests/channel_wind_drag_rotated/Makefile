PROJECT = channel

input: clean
	gmsh -3 -o channel.msh src/channel.geo
	../../bin/gmsh2triangle channel.msh
	../../bin/transform_mesh '(cos(1.)*x-sin(1.0)*y,sin(1.)*x+cos(1.)*y,z)' channel
        

clean:
	rm -f $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node *.bak
	rm -f *.msh
