FLMLMODEL1 = top_hat_cg_explicit_pe
FLMLMODEL2 = top_hat_cg_su

default: input

input:

clean: clean-run-debug
clean-run:
	rm -f $(FLMLMODEL1)_?*.pvtu $(FLMLMODEL1)_?*.vtu
	rm -f $(FLMLMODEL1).detectors $(FLMLMODEL1).stat tmpf*
	rm -f $(FLMLMODEL1)_?*_checkpoint.flml
	rm -f $(FLMLMODEL1)_?*_checkpoint_?*.ele $(FLMLMODEL1)_?*_checkpoint_?*.bound $(FLMLMODEL1)_?*_checkpoint_?*.halo $(FLMLMODEL1)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL1)_?*_?*_checkpoint.pvtu $(FLMLMODEL1)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL1)_?*_checkpoint.ele $(FLMLMODEL1)_?*_checkpoint.bound $(FLMLMODEL1)_?*_checkpoint.node
	rm -f $(FLMLMODEL1)_?*_checkpoint.vtu
	rm -f $(FLMLMODEL2)_?*.pvtu $(FLMLMODEL2)_?*.vtu
	rm -f $(FLMLMODEL2).detectors $(FLMLMODEL2).stat tmpf*
	rm -f $(FLMLMODEL2)_?*_checkpoint.flml
	rm -f $(FLMLMODEL2)_?*_checkpoint_?*.ele $(FLMLMODEL2)_?*_checkpoint_?*.bound $(FLMLMODEL2)_?*_checkpoint_?*.halo $(FLMLMODEL2)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL2)_?*_?*_checkpoint.pvtu $(FLMLMODEL2)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL2)_?*_checkpoint.ele $(FLMLMODEL2)_?*_checkpoint.bound $(FLMLMODEL2)_?*_checkpoint.node
	rm -f $(FLMLMODEL2)_?*_checkpoint.vtu
	rm -f matrixdump matrixdump.info
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.bound adapted_mesh.ele adapted_mesh.node
	rm -f adapted_mesh_?*.bound adapted_mesh_?*.ele adapted_mesh_?*.node
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
