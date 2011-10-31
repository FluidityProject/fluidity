input: clean
	../../bin/gmsh2triangle --2d src/channel.msh
	../../bin/fldecomp -n 2  channel

clean:
	rm -fr  *.ele *.edge *.node *.poly *.halo *.vtu *.pvtu fluidity.*  *.s *.d.1 *.stat tracer_*
