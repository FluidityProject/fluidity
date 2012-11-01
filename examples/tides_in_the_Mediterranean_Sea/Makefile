include envcheck.mk

preprocess: envcheck
ifeq ("x","x${NPROCS}")
	@echo "**********ERROR: Please define NPROCS to set the number of processors to use"
	@exit 1
endif
	@echo **********Converting the gmsh mesh to triangle format:
	${BINPREFIX}gmsh2triangle med.msh
	@echo **********Modifying mesh nodal file using sed:
	sed '1s/3/2/2' < med.node > med2.node
	@echo **********Replacing mesh with modified version:
	mv med2.node med.node
ifneq ($(NPROCS),64)
	@echo **********WARNING: This simulation is best run on 64 processors
endif
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
else 
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	${BINPREFIX}fldecomp -n $(NPROCS) -s med 
endif

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
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}fluidity -v2 -l med.flml
endif

postprocess: envcheck
	@echo **********Reprojecting final output to long-lat and extracting RMS differences:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python ./reproject-to-loglat.py
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python ./Med-tides-probe.py


clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf  *gmt* *.ele *.node *.edge Med-GEBCO-5m-gauges-fes2004-O1-102 *med_* *.halo *.pvtu *.vtu *.stat *.log-* *.err-* *check* *Check* *matrixdump*

input: clean
	$(MAKE) preprocess NPROCS=64

