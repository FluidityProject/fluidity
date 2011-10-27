input: clean
	../../bin/create_aligned_mesh box 50 50 -- -50 10 10 10
	../../bin/fldecomp -n 2 box

clean:
	# serial
	rm -f  *.ele *.edge *.face *.msh *.node *.poly *.vtu *.s *.d.1 *.stat *.detectors*
	rm -f matrix*
	rm -f fluidity.*
	# parallel
	rm -rf biology_conservation_parallel_*
	rm -rf *.halo
