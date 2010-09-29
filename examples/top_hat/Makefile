preprocess: 
	@echo **********Generating a 1D mesh for the [0,3] interval:
	../../scripts/interval --dx=0.025 --reverse 0 3 line

run:
	@echo **********Calling fluidity with verbose log output enabled:
	../../bin/fluidity -v2 -l top_hat_cg_adapt.flml

postprocess:
	@echo **********No preprocessing needed

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf line.* *.stat *.node *.ele *.bound top_hat*vtu *.log-0 *.err-0 matrixdump*
