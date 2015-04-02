include envcheck.mk
NPROCS=1

preprocess: envcheck
	install -d "simulation"; \
	cd simulation; \
	ln -s ../backward_facing_step_2d.flml .; \
	echo **********Generating mesh; \
	gmsh -2 -o step2d-kim-med.msh ../src/step2d-kim-med.geo; \
	cd ..

run: envcheck
	@echo **********Running simulation:
	cd simulation &&  ${BINPREFIX}fluidity -v2 -l backward_facing_step_2d.flml &&  cd ..

postprocess: envcheck
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	cd simulation; \
	ln -s ../postprocessor_2d.py .; \
	./postprocessor_2d.py; \
	cd ..

input: clean
	$(MAKE) preprocess 

clean:
	@echo **********Cleaning the output from previous fluidity runs
	rm -rf simulation *.pdf *flredecomp*

