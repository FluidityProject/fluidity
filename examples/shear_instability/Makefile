include envcheck.mk

ifndef NPROCS
NPROCS=1
endif	

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 -v 0 tank_fixed_mesh.geo

	@echo **********Generating the input files:
	export SIMULATION_NAME=shear_kelvin_helmholtz_low_res; \
	cat shear.flml.in | EXPANSION_COEFF=0.0 MIN_LENGTH=0.3 envsubst > $$SIMULATION_NAME.flml
	export SIMULATION_NAME=shear_holmboe_low_res; \
	cat shear.flml.in | EXPANSION_COEFF=1.0e-3 MIN_LENGTH=0.3 envsubst > $$SIMULATION_NAME.flml
	export SIMULATION_NAME=shear_kelvin_helmholtz_high_res; \
	cat shear.flml.in | EXPANSION_COEFF=0.0 MIN_LENGTH=0.1 envsubst > $$SIMULATION_NAME.flml
	export SIMULATION_NAME=shear_holmboe_high_res; \
	cat shear.flml.in | EXPANSION_COEFF=1.0e-3 MIN_LENGTH=0.1 envsubst > $$SIMULATION_NAME.flml
	export SIMULATION_NAME=shear_kelvin_helmholtz_fixed_mesh; \
	cat shear_fixed.flml.in | EXPANSION_COEFF=0.0 envsubst > $$SIMULATION_NAME.flml
	export SIMULATION_NAME=shear_holmboe_fixed_mesh; \
	cat shear_fixed.flml.in | EXPANSION_COEFF=1.0e-3 envsubst > $$SIMULATION_NAME.flml

	@echo **********Running periodise
	${BINPREFIX}periodise shear_kelvin_helmholtz_low_res.flml
	${BINPREFIX}periodise shear_holmboe_low_res.flml	
	${BINPREFIX}periodise shear_kelvin_helmholtz_high_res.flml
	${BINPREFIX}periodise shear_holmboe_high_res.flml
	${BINPREFIX}periodise shear_kelvin_helmholtz_fixed_mesh.flml
	${BINPREFIX}periodise shear_holmboe_fixed_mesh.flml
ifneq ($(NPROCS),1)
	mpiexec -n ${NPROCS} ${BINPREFIX}flredecomp -i 1 -o ${NPROCS} shear_kelvin_helmholtz_fixed_mesh_periodised shear_kelvin_helmholtz_fixed_mesh_periodised_redecomped
	mpiexec -n ${NPROCS} ${BINPREFIX}flredecomp -i 1 -o ${NPROCS} shear_holm_boe_fixed_mesh_periodised shear_holmboe_fixed_mesh_periodised_redecomped
else
	@echo **********WARNING: The fixed mesh run may take a long time to run in serial
	@echo **********Set NPROCS to run in parallel.
endif

run: envcheck preprocess
ifeq ($(NPROCS),1)
	@echo **********WARNING: The fixed mesh run may take a long time to run in serial
	@echo ********** Set NPROCS to run in parallel.
endif
	@echo **********Calling fluidity with least verbose log output enabled:
	${BINPREFIX}fluidity -v0 -l shear_kelvin_helmholtz_low_res_periodised.flml
	${BINPREFIX}fluidity -v0 -l shear_holmboe_low_res_periodised.flml
	${BINPREFIX}fluidity -v0 -l shear_kelvin_helmholtz_high_res_periodised.flml
	${BINPREFIX}fluidity -v0 -l shear_holmboe_high_res_periodised.flml
ifeq ($(NPROCS),1)
	mpiexec -n ${NPROCS} ${BINPREFIX}fluidity -v0 -l shear_kelvin_helmholtz_fixed_mesh_periodised_redecomped.flml
	mpiexec -n ${NPROCS} ${BINPREFIX}fluidity -v0 -l shear_holmboe_fixed_mesh_periodised_redecomped.flml
else
	${BINPREFIX}fluidity -v0 -l shear_kelvin_helmholtz_fixed_mesh_periodised.flml
	${BINPREFIX}fluidity -v0 -l shear_holmboe_fixed_mesh_periodised.flml
endif

postprocess: envcheck

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.msh *.node *.ele *.edge shear*.vtu *checkpoint* *convergence* *.log-* *.err-* matrixdump* *periodised* shear_fixed_mesh_[0-9]* *.flml

