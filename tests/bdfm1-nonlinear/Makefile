input: clean
	gmsh -3 channel.geo
	../../tools/gmsh2triangle.py channel.msh
	gmsh -3 channel_unstr.geo
	../../tools/gmsh2triangle.py channel_unstr.msh

clean:
	rm -f *.d.* *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f channel*.ele channel*.face channel*.node channel*.edge *.stat *.halo
	rm -f *.msh
