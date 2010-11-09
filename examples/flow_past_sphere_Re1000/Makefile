# Default number of processes:
NPROCS=1
# directory locations:
BIN=../../bin/
SCRIPTS=../../scripts/

preprocess:
	@echo **********Converting the gmsh mesh to triangle format:
	$(SCRIPTS)gmsh2triangle unit_sphere.msh
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
else
# exit if no fltools and print warning
	@if [ ! -e $(BIN)fldecomp ] ; then echo "ERROR: You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	$(BIN)fldecomp -n $(NPROCS) -f unit_sphere
endif


run:
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	$(BIN)fluidity -v2 -l flow_past_sphere_Re1000.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) $(BIN)fluidity -v2 -l flow_past_sphere_Re1000.flml
endif

postprocess:
	@echo **********No postprocessing needed

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *vtu *.stat *.log-? *.err-? matrixdump* 
	rm -rf flow_past_sphere_Re1000_? flow_past_sphere_Re1000_??
