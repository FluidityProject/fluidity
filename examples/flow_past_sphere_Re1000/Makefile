include envcheck.mk

preprocess: envcheck
	@echo **********Converting the gmsh mesh to triangle format:

run: envcheck
ifneq ($(NPROCS),8)
	@echo **********WARNING: This simulation is best run on 8 processors
endif
ifeq ($(NPROCS),1)
	@echo **********WARNING: This is a large simulation and will take a very long time in serial. Find a handy supercomputer.
	@echo **********Calling fluidity in serial with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l flow_past_sphere_Re1000.flml
else
	@echo **********Calling flredecomp in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}flredecomp -i 1 -o $(NPROCS) -v -l flow_past_sphere_Re1000 flow_past_sphere_Re1000_flredecomp
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) ${BINPREFIX}fluidity -v2 -l flow_past_sphere_Re1000_flredecomp.flml
endif

postprocess:
	@echo **********Calling the data extraction and plotting script
	python ./plot_data.py

input: clean
	$(MAKE) preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *vtu *.stat *.log-* *.err-* matrixdump* 
	rm -rf flow_past_sphere_Re1000_? flow_past_sphere_Re1000_?? 
	rm -rf first_timestep_adapted_mesh*
	rm -f Sphere_drag.pdf
	rm -rf *flredecomp*
