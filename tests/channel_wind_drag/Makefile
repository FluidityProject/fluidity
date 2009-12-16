PROJECT = channel

input: clean
	gmsh -3 -o channel.msh src/channel.geo
	../../scripts/gmsh2triangle channel.msh

clean:
	rm -f $(PROJECT) $(PROJECT)0 $(PROJECT).fld
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node 
	rm -f *.msh
