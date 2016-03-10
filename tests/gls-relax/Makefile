input: clean
	# note the -- -50 - we want the z axis with 0 at surface and -50m at depth
	../../bin/create_aligned_mesh gls-MixedLayer 50 50 -- -50 3 3 5

clean: 
	rm -f  *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *log *err
	rm -f *.msh fluidity.log_*
