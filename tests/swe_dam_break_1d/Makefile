input: clean preprocess

preprocess:
	../../bin/interval --dx=5 -- 0 2000 line

run:
	../../bin/fluidity -v2 -l swe_dam_break_1d.flml
	
clean:
	rm -f *.node *.ele *.bound
	rm -f *.vtu *.stat
	rm -f fluidity.err-0 fluidity.log-0
	rm -f matrixdump*
