preprocess: 
	@echo **********Generating a 1D mesh for the [0,3] interval:
	../../bin/interval --dx=0.025 --reverse 0 3 line

run:
	@echo **********Running the Continuous Galerkin version of this example:
	../../bin/fluidity -v2 -l top_hat_cg.flml
	@echo **********Running the Discontinuous Galerkin version of this example:
	../../bin/fluidity -v2 -l top_hat_dg.flml
	@echo **********Running the Control Volumes version of this example:
	../../bin/fluidity -v2 -l top_hat_cv.flml

postprocess:
	@echo **********Calling visualise.py
	python visualise.py top_hat_cg_100.vtu
	python visualise.py top_hat_cv_100.vtu
	python visualise.py top_hat_dg_100.vtu

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf line.* *.stat *.node *.ele *.bound top_hat*vtu *.log-0 *.err-0 matrixdump*
