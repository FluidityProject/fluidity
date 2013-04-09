SIMULATION=sphericalShell
MESH_FILE_PATTERN=sphericalShell
FLUIDITY=../../


BIN=$(FLUIDITY)bin
SCHEMAS=$(FLUIDITY)schemas
#Export MESH_FILE_PATTERN, in order to suply variable to
#mesh/Makefile
export SIMULATION
#export MESH_FILE_PATTERN
#export FLUIDITY
#export BIN
PYTHONPATH := $(FLUIDITY)python/:$(BIN)
FLML_FILE=$(SIMULATION).flml
STAT_FILE=$(SIMULATION).stat
VTU_FILE_PATTERN=$(SIMULATION)_..vtu
PVTU_FILE_PATTERN=$(SIMULATION)_..pvtu

all: run

run: input
ifeq ($(NPROCS),1)
	time $(BIN)/fluidity -v5 -l $(FLML_FILE) &
else
	time mpiexec -n $(NPROCS) $(BIN)/fluidity -v3 -l $(FLML_FILE) &
endif

gdb_run: input
ifeq ($(NPROCS),1)
	gdb --args $(BIN)/fluidity -v5 -l $(FLML_FILE) &
else
	#mpirun -n $(NPROCS) xterm -e gdb --batch -x gdbCommands.txt --args $(BIN)/fluidity -v5 -l $(FLML_FILE) &
	mpirun -n $(NPROCS) xterm -e gdb --args $(BIN)/fluidity -v5 -l $(FLML_FILE)
endif

input: $(FLML_FILE) mesh
ifndef NPROCS
	@echo "**********ERROR: NPROCS not defined."
	@false
endif
ifeq ($(NPROCS),1)
	@echo "**********WARNING: NPROCS is equal to 1"
	@echo "**********WARNING: Serial run."
	@echo "**********WARNING: This is a large simulation and will take"
	@echo "                   a very long time in serial."
endif
	$(MAKE) -e -C mesh/ mesh
ifeq ($(NPROCS),1)
	@echo "**********WARNING: Not decomposing mesh."
else
# exit if no fltools and print warning.
	@if [ ! -e $(BIN)/flredecomp ] ; then echo "**********ERROR: You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo "**********Decomposing the mesh into ${NPROCS} parts for parallel run:"
	$(BIN)/fldecomp -n ${NPROCS} -s mesh/$(MESH_FILE_PATTERN)
endif

clean:clean-light
	@rm -f *.s *.d.1 *.stat where_am_i
	$(MAKE) -e -C mesh/ $@

clean-light:clean-results clean-logs clean-stat clean-matrixdump clean-checkpoints clean-pyc

clean-results:
	@echo "**********Cleaning results ..."
	@rm -rf $(SIMULATION)_[0-9]* *.vtu *.pvtu

clean-logs:
	@echo "**********Cleaning log-files ..."
	@rm -f  *.err-* *.log-*

clean-grid:
	@echo "**********Cleaning mesh-files ..."
	@$(MAKE) -C mesh/ clean

clean-stat:
	@echo "**********Cleaning stat and detector files ..."
	@rm -f  *.stat

clean-matrixdump:
	@echo "**********Cleaning matrixdump files ..."
	@rm -f matrixdump*

clean-checkpoints:
	@echo "**********Cleaning checkpoints ..."
	@rm -rf *checkpoint*

clean-pyc:
	@echo "**********Cleaning python compiled files ..."
	@rm -rf *.pyc

config:
	diamond -s $(SCHEMAS)/fluidity_options.rng  $(FLML_FILE) &

viz:
	paraview --data=$(VTU_FILE_PATTERN) 2>/dev/null &

vizpara:
	paraview --data=$(PVTU_FILE_PATTERN) 2>/dev/null &

vizstat:
	$(BIN)/statplot $(STAT_FILE) &

.PHONY: mesh
