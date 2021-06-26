include envcheck.mk
NPROCS=1

preprocess: envcheck
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make preprocess TYPE=type, where type is reference or kepsilon"; false; fi;
	echo **********Installing directory for $(TYPE) run:
	install -d $(TYPE); \
	cd $(TYPE); \
	ln -s ../backward_facing_step_2d_kim_$(TYPE).flml .; \
	echo **********Generating mesh for $(TYPE) simulation.; \
	gmsh -2 -o step2d-kim-med.msh ../src/step2d-kim-med.geo; \
	cd ..

run: envcheck
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make run TYPE=type, where type is reference or kepsilon"; false; fi;
	@echo **********Running the $(TYPE) simulation:
	cd $(TYPE) &&  ${BINPREFIX}fluidity -v2 -l backward_facing_step_2d_kim_$(TYPE).flml &&  cd ..

postprocess: envcheck
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make postprocess TYPE=type, where type is reference or kepsilon"; false; fi;
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	cd $(TYPE); \
	ln -s ../postprocessor_2d.py .; \
	./postprocessor_2d.py $(TYPE); \
	cd ..

input: clean
	$(MAKE) preprocess TYPE=reference
	$(MAKE) preprocess TYPE=kepsilon

clean:
	@echo **********Cleaning the output from previous fluidity runs
	rm -rf kepsilon reference *.pdf *flredecomp*

