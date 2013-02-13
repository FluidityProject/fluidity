include envcheck.mk

preprocess: envcheck
	@echo **********Generating a 1D mesh for the [0,3] interval:
	${BINPREFIX}interval --dx=0.025 --reverse 0 3 line
	${BINPREFIX}interval --dx=0.025 --reverse 0 3 line

run: envcheck
	@echo **********Running the Continuous Galerkin version of this example:
	${BINPREFIX}fluidity -v2 -l top_hat_cg.flml
	@echo **********Running the Discontinuous Galerkin version of this example:
	${BINPREFIX}fluidity -v2 -l top_hat_dg.flml
	@echo **********Running the Control Volumes version of this example:
	${BINPREFIX}fluidity -v2 -l top_hat_cv.flml

postprocess: envcheck
	@echo **********Calling visualise.py
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python visualise.py top_hat_cg_100.vtu
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python visualise.py top_hat_cv_100.vtu
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python visualise.py top_hat_dg_100.vtu

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf line.* *.stat *.node *.ele *.bound top_hat*vtu *.log-0 *.err-0 matrixdump*
