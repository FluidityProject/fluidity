FLMLMODEL1 = top_hat_cg_explicit_pe
FLMLMODEL2 = top_hat_cg_su

default: input

input:
	../../bin/interval --dx=0.025 0 1 line

clean: clean-run-debug
	rm -f *.log
clean-run:
	rm -f $(FLMLMODEL1)_?*.pvtu $(FLMLMODEL1)_?*.vtu
	rm -f $(FLMLMODEL1).detectors $(FLMLMODEL1).stat tmpf*
	rm -f $(FLMLMODEL1)_?*_checkpoint.flml
	rm -f $(FLMLMODEL1)_?*_checkpoint_?*.msh
	rm -f $(FLMLMODEL1)_?*_?*_checkpoint.pvtu $(FLMLMODEL1)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL1)_?*_checkpoint.msh
	rm -f $(FLMLMODEL1)_?*_checkpoint.vtu
	rm -f $(FLMLMODEL2)_?*.pvtu $(FLMLMODEL2)_?*.vtu
	rm -f $(FLMLMODEL2).detectors $(FLMLMODEL2).stat tmpf*
	rm -f $(FLMLMODEL2)_?*_checkpoint.flml
	rm -f $(FLMLMODEL2)_?*_checkpoint_?*.msh
	rm -f $(FLMLMODEL2)_?*_?*_checkpoint.pvtu $(FLMLMODEL2)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL2)_?*_checkpoint.msh
	rm -f $(FLMLMODEL2)_?*_checkpoint.vtu
	rm -f matrixdump matrixdump.info
clean-run-debug: clean-run
	rm -f line.msh
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.msh
	rm -f adapted_mesh_?*.msh
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
