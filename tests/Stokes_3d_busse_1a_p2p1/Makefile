input: clean
	gmsh -3 -optimize Busse_C1.msh src/Busse_C1.geo 
	../../scripts/gmsh2triangle Busse_C1.msh

clean:
	rm -f *.ele *.edge *.node *.vtu *.stat *.msh *.detectors *.dat fluidity.* *_state

