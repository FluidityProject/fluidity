FLMLMODEL = square-convection
MESH = square

default: input

input: clean
	tar -xvjf src/$(MESH).tar.bz2
	triangle2gmsh $(MESH)

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).ele $(MESH).edge $(MESH).node $(MESH).msh
	rm -f $(MESH)_?*.msh $(MESH)_?*.halo decomp.log decomp_err.log
clean-run:
	rm -rf $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*.vtu $(FLMLMODEL)_?*
	rm -f $(FLMLMODEL).detectors $(FLMLMODEL).detectors.dat $(FLMLMODEL).stat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f $(FLMLMODEL)_?*_checkpoint.flml
	rm -f $(FLMLMODEL)_?*_checkpoint_?*.msh $(FLMLMODEL)_?*_checkpoint_?*.halo
	rm -f $(FLMLMODEL)_?*_?*_checkpoint.pvtu $(FLMLMODEL)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL)_?*_checkpoint.msh $(FLMLMODEL)_?*_checkpoint.halo
	rm -f $(FLMLMODEL)_?*_checkpoint.vtu
	rm -rf *flredecomp*
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.msh adapted_mesh.halo
	rm -f adapted_mesh_?*.msh adapted_mesh_?*.halo
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
