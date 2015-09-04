input: clean preprocess

preprocess:
	gmsh -2 src/domain.geo -o src/domain.msh

run:
	../../bin/fluidity swe_tidal_flow_regular_bed.flml
	
clean:
	rm -f *.node *.ele *.bound
	rm -f *.vtu *.stat
	rm -f fluidity.err-0 fluidity.log-0
	rm -f matrixdump*
