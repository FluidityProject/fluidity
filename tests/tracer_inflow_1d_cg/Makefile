input: clean
	# ../../bin/interval --dx 0.5 0 10 interval
	gmsh -1 interval.geo

clean:
	rm -f  *.ele *.bound *.node *.poly *.vtu *.s *.d.1 *.stat \
	matrixdump matrixdump.info
