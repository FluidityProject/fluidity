include envcheck.mk

# Default number of processes:
NPROCS=1

preprocess: envcheck
ifeq ("x","x${NPROCS}")
	@echo "**********ERROR: Please define NPROCS to set the number of processors to use"
	@exit 1
endif
	@echo **********Converting the gmsh mesh to triangle format:
	${BINPREFIX}gmsh2triangle unit_sphere.msh
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
else
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	${BINPREFIX}fldecomp -n $(NPROCS) -f unit_sphere
endif

run: envcheck
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l flow_past_sphere_Re100.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}fluidity -v2 -l flow_past_sphere_Re100.flml
endif

postprocess:
	@echo **********Calling the data extraction and plotting script
	python ./plot_data.py

input: clean
	$(MAKE) preprocess NPROCS=8

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *vtu *.stat *.log-? *.err-? matrixdump* 
	rm -rf flow_past_sphere_Re100_? flow_past_sphere_Re100_?? 
	rm -f Sphere_drag.pdf
