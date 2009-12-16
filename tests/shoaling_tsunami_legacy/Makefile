input: clean
	gmsh -3 src/shoaling_tsunami.geo -o shoaling_tsunami.msh
	../../scripts/gmsh2triangle shoaling_tsunami.msh
	src/bathy.py
	mv output.node shoaling_tsunami.node

clean:
	rm -rf *.vtu *.stat *.detectors *.msh *.node *.ele *.face *log* *err* matrixdump*
