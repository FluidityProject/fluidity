include envcheck.mk

preprocess: envcheck

run: envcheck
ifndef ($(NPROCS))
	@echo **********ERROR: Please define NPROCS to set the number of processors to use
	@exit 1
endif
ifneq ($(NPROCS),64)
	@echo **********WARNING: This simulation is best run on 64 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l med.flml
else
	@echo **********Calling flredecomp in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}flredecomp -i 1 -o $(NPROCS) -v -l med med_flredecomp
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}fluidity -v2 -l med_flredecomp.flml
endif

postprocess: envcheck
	@echo **********Reprojecting final output to long-lat and extracting RMS differences:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python ./reproject-to-loglat.py
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python ./Med-tides-probe.py


clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf  *gmt* *.ele *.node *.edge Med-GEBCO-5m-gauges-fes2004-O1-102 *med_* *.halo *.pvtu *.vtu *.stat *.log-* *.err-* *checkp* *Check* *matrixdump*
	rm -rf *flredecomp*

input: clean
	$(MAKE) preprocess

