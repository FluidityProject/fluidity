include envcheck.mk

NPROCS?=1

preprocess: envcheck
	@echo **********Copying circle.geo to the current directory:
	cp src/circle.geo .
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 circle.geo 
	@echo **********Converting the gmsh mesh to triangle format:

run: envcheck
ifneq ($(NPROCS),64)
	@echo **********WARNING: This simulation is best run on 64 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l restratification_after_oodc.flml
else
	@echo **********Calling flredecomp in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}flredecomp -i 1 -o $(NPROCS) -v -l restratification_after_oodc restratification_after_oodc_flredecomp
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}fluidity -v2 -l restratification_after_oodc_flredecomp.flml
endif

postprocess:
	@echo **********No postprocessing needed

input: clean
	$(MAKE) preprocess NPROCS=64

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.geo *.msh *.ele *.edge *.face *.node *.halo *.poly *.vtu *.pvtu *.stat *.log-* *.err-* */matrixdump*
	rm -rf restratification_after_oodc_? restratification_after_oodc_??
	rm -rf *flredecomp*

realclean: clean
	@echo **********Cleaning the output and checkpoints from previous fluidity runs:
	rm -rf *checkpoint*
