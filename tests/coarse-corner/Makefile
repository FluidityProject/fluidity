FLMLMODEL = coarse-corner
MESH = CoarseCorner

default: input

input: src/$(MESH).geo clean
	gmsh -3 src/$(MESH).geo -o $(MESH).msh

clean: clean-mesh clean-run-debug clean-checkpoint
clean-mesh:
	rm -f $(MESH).msh
	rm -f $(MESH)_?*.msh $(MESH)_?*.halo
	rm -f decomp.log decomp_err.log

clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?/$(FLMLMODEL)_?*.vtu
	rm -rf $(FLMLMODEL)_?
	rm -f $(FLMLMODEL).detectors $(FLMLMODEL).stat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f $(FLMLMODEL)_?*_checkpoint.flml
	rm -f $(FLMLMODEL)_?*_checkpoint_?*.msh $(FLMLMODEL)_?*_checkpoint_?*.halo
	rm -f $(FLMLMODEL)_?*_?*_checkpoint.pvtu $(FLMLMODEL)_?*_?*_checkpoint/$(FLMLMODEL)_?*_?*_checkpoint.vtu
	rm -rf $(FLMLMODEL)_?*_?*_checkpoint
	rm -f $(FLMLMODEL)_?*_checkpoint.msh $(FLMLMODEL)_?*_checkpoint.halo
	rm -f $(FLMLMODEL)_?*_checkpoint.vtu
	rm -rf *flredecomp*
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.msh adapted_mesh.halo
	rm -f adapted_mesh_?*.msh adapted_mesh_?*.halo
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
clean-checkpoint:
	make -C checkpoint clean
