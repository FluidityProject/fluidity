include envcheck.mk

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 3d:
	gmsh -3 -optimize -o step3d.msh src/step3d.geo

run: envcheck
	@if [ "x" == "x${NPROCS}" ]; then echo "**********ERROR: Please define NPROCS to set the number of processors to use"; false; fi
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l backward_facing_step_3d.flml
else
	@echo **********Calling flredecomp in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}flredecomp -i 1 -o $(NPROCS) -v -l backward_facing_step_3d backward_facing_step_3d_flredecomp
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}fluidity -v2 -l backward_facing_step_3d_flredecomp.flml
endif

postprocess: envcheck
	@echo **********Calling the velocity data extraction and plotting python script:
	@mkdir -p numpy_data 
	./postprocessor_3d.py

input: clean
	$(MAKE) preprocess NPROCS=8

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.stat *.msh backward_facing_step_3d_*.pvtu backward_facing_step_3d_* *.log-* *.err-* matrixdump* *.halo
	rm -rf *flredecomp*
