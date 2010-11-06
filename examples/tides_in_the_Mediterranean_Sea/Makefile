# Default number of processes (64 processors):
NPROCS=64
# directory locations:
BIN=../../bin/
SCRIPTS=../../scripts/

preprocess: 
	@echo **********Converting the gmsh mesh to triangle format:
	$(SCRIPTS)gmsh2triangle med.msh
	@echo **********Modifying mesh nodal file using sed:
	sed < med.node 's/[0-9]*/2/2' > med2.node
	@echo **********Replacing mesh with modified version:
	mv med2.node med.node
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
else
	@if [ ! -e ../../bin/fldecomp ] ; then echo "ERROR: You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	../../bin/fldecomp -n $(NPROCS) -s med 
endif

run:
ifeq ($(NPROCS),1)
	@echo **********Calling fluidity in serial with verbose log output enabled:
	$(BIN)fluidity -v2 -l med.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) $(BIN)fluidity -v2 -l med.flml
endif

postprocess:
	@echo **********Reprojecting final output to long-lat and extracting RMS differences:
	python ./reproject-to-loglat.py
	python ./Med-tides-probe.py


clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf  *gmt* *.ele *.node *.edge Med-GEBCO-5m-gauges-fes2004-O1-102 *med_* *.halo *.pvtu *.vtu *.stat *.log-* *.err-* *check* *Check* *matrixdump*

input: clean preprocess
