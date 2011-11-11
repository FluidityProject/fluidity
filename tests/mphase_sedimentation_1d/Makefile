preprocess: 
	@echo **********Generating a 1D mesh for the [0,1] interval:
	../../bin/interval --dx=0.025 --reverse 0 1 line

run:
	@echo **********Running simulation using 1D mesh
	../../bin/fluidity -v2 -l mphase_sedimentation_1d.flml

input: clean preprocess

clean:
	rm -f *.stat *.steady_state*
	rm -f *.d.* *.vtu
	rm -f *.msh
	rm -f *.ele *.edge *.node *.poly *.bound
	rm -f matrixdump* *log* *err*

