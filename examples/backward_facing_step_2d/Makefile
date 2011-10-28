BIN=../../../bin
preprocess:
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make preprocess TYPE=type, where type is reference or kepsilon"; false; fi;
	echo **********Installing directory for $(TYPE) run:
	install -d $(TYPE); \
	cd $(TYPE); \
	ln -s ../backward_facing_step_2d_kim_$(TYPE).flml .; \
	echo **********Generating mesh for $(TYPE) simulation.; \
	gmsh -2 -o step2d-kim-med.msh ../src/step2d-kim-med.geo; \
	$(BIN)/gmsh2triangle -2 step2d-kim-med.msh; \
	cd ..

run:
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make run TYPE=type, where type is reference or kepsilon"; false; fi;
	@echo **********Running the $(TYPE) simulation:
	cd $(TYPE) &&  $(BIN)/fluidity -v2 -l backward_facing_step_2d_kim_$(TYPE).flml &&  cd ..

postprocess:
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make postprocess TYPE=type, where type is reference or kepsilon"; false; fi;
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	cd $(TYPE); \
	ln -s ../postprocessor_2d.py .; \
	./postprocessor_2d.py $(TYPE); \
	cd ..

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	rm -r kepsilon reference *.pdf

