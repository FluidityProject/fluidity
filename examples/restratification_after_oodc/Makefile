# Default number of processes (32 processors):
NPROCS=32
# directory locations:
BIN=../../bin/
SCRIPTS=../../scripts/

preprocess:
	@echo **********Copying circle.geo to the current directory:
	cp src/circle.geo .
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 circle.geo 
	@echo **********Converting the gmsh mesh to triangle format:
	$(SCRIPTS)gmsh2triangle --2d circle.msh; 
ifeq ($(NPROCS),1)
	@echo **********Serial run: not decomposing mesh
else
# exit if no fltools and print warning
	@if [ ! -e $(BIN)fldecomp ] ; then echo "ERROR: You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(NPROCS) parts for parallel run:
	$(BIN)fldecomp -n $(NPROCS) -f circle
endif

run:
ifeq ($(NPROCS),1)
	@echo **********Calling fluidity in serial with verbose log output enabled:
	$(BIN)fluidity -v2 -l restratification_after_oodc.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n $(NPROCS) $(BIN)fluidity -v2 -l restratification_after_oodc.flml
endif
        
postprocess:
	@echo **********No postprocessing needed

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *.vtu *.stat *.log-0 *.err-0 */matrixdump*
