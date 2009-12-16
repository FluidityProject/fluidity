PROJECT = channel

input: clean
	gmsh -3 -o channel.msh src/channel.geo
	../../scripts/gmsh2triangle channel.msh
	../../scripts/transform_mesh '(x,0.0 if y==0.0 else exp(-10.0*log(0.5)/1.e5*y)/exp(-10.0*log(0.5))*1.e5,z)' channel

clean:
	rm -f $(PROJECT) $(PROJECT)0 $(PROJECT).fld
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node *.msh
	rm -f matrixdump*
