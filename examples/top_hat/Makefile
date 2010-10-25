preprocess: 
	@echo **********Generating a 1D mesh for the [0,3] interval:
	../../scripts/interval --dx=0.025 --reverse 0 3 line

run:
	@echo **********Running the Continuous Galerkin version of this example:
	../../bin/fluidity -v2 -l top_hat_cg.flml
	@echo **********Running the Discontinuous Galerkin version of this example:
	../../bin/fluidity -v2 -l top_hat_dg.flml
	@echo **********Running the Control Volumes version of this example:
	../../bin/fluidity -v2 -l top_hat_cv.flml

postprocess:
	@echo **********No postprocessing needed

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf line.* *.stat *.node *.ele *.bound top_hat*vtu *.log-0 *.err-0 matrixdump*
