input: clean
	gmsh -3 src/cube-1.geo 
	../../bin/gmsh2triangle src/cube-1.msh

clean:
	rm -f src/cube-1.msh *.ele *.face *.node *.vtu *.s *.d.1 *.stat *.log-0 *.err-0
