input: clean
	../../scripts/gmsh2triangle -2 src/circle-serial.msh
	../../scripts/gmsh2triangle -2 src/circle-parallel.msh
	../../bin/fldecomp -n 3  circle-parallel

clean:
	rm -f *.halo *.ele *.edge *.face *.node *.poly *.vtu *.pvtu *.s *.d.* *.stat *.s *.log-0 *.err-0
