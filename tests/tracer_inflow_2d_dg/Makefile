input: clean
	../../scripts/gmsh2triangle --2d src/channel.msh
	../../bin/fldecomp -n 2  channel

clean:
	rm -f  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.*  *.s *.d.1 *.stat
