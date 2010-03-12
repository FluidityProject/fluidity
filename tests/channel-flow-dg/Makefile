PROJECT = channel

input: clean
	../../scripts/gmsh2triangle --2d ./src/channel.msh

clean:
	rm -f $(PROJECT) $(PROJECT)0 $(PROJECT).fld
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node 
	rm -f *.msh
