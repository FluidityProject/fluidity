input: clean preprocess

preprocess:
	gmsh -2 src/domain.geo -o src/domain.msh

run:
	../../bin/fluidity -v2 -l swe_dam_break_2d.flml
	
clean:
	rm -f *.node *.ele *.bound
	rm -f *.vtu *.stat
	rm -f fluidity.err-0 fluidity.log-0
	rm -f matrixdump*
