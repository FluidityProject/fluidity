BIN=../../../bin
preprocess:
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make preprocess TYPE=type, where type is reference or kepsilon"; false; fi;
	echo **********Installing directory for $(TYPE) run:
	install -d $(TYPE)-$(MESH); \
	cd $(TYPE)-$(MESH); \
	ln -s ../flmls/backward_facing_step_2d-$(TYPE)-$(MESH).flml .; \
	echo **********Generating mesh for $(TYPE)-$(MESH) simulation.; \
	gmsh -2 -optimize -o step2d-kim-$(MESH).msh ../src/step2d-kim-$(MESH).geo; \
	$(BIN)/gmsh2triangle -2 step2d-kim-$(MESH).msh; \
	cd ..

run:
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make run TYPE=type, where type is reference or kepsilon"; false; fi;
	@echo **********Running the $(TYPE)-$(MESH) simulation:
	cd $(TYPE)-$(MESH) &&  $(BIN)/fluidity -v2 -l backward_facing_step_2d-$(TYPE)-$(MESH).flml &&  cd ..

postprocess:
	@if [ "x" == "x${TYPE}" ]; then echo "**********ERROR: Commandline should be make postprocess TYPE=type, where type is reference or kepsilon"; false; fi;
	@echo **********Calling the velocity data extraction and reattachment length calculation python script:
	cd $(TYPE)-$(MESH); \
	ln -s ../postprocessor_2d.py .; \
	./postprocessor_2d.py $(TYPE) $(MESH); \
	cd ..

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs
	rm -r kepsilon reference *.pdf

space:
	@$(foreach N,1 2 3 4 5 6 7 8 9, \
		rm */backward_facing_step_2d*_*$(N).vtu; \
	)
