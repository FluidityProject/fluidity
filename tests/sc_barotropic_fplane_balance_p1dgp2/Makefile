FLMLMODEL1 = sc_barotropic_fplane_balance_p1dgp2_0
FLMLMODEL2 = sc_barotropic_fplane_balance_p1dgp2_1
MESH = square

default: input

input:
	gmsh -2 src/$(MESH).geo
	cp src/$(MESH).msh .

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f src/$(MESH).msh
	rm -f $(MESH).edge $(MESH).ele $(MESH).node
	rm -f $(MESH)_?*.edge  $(MESH)_?*.ele $(MESH)_?*.halo $(MESH)_?*.node decomp.log decomp_err.log
clean-run:
	rm -f matrixdump matrixdump.info
	rm -f $(FLMLMODEL1)_?*.pvtu $(FLMLMODEL1)_?*.vtu
	rm -f $(FLMLMODEL1).detectors $(FLMLMODEL1).stat tmpf*
	rm -f $(FLMLMODEL1).log
	rm -f $(FLMLMODEL1)_?*_checkpoint.flml
	rm -f $(FLMLMODEL1)_?*_checkpoint_?*.ele $(FLMLMODEL1)_?*_checkpoint_?*.face $(FLMLMODEL1)_?*_checkpoint_?*.halo $(FLMLMODEL1)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL1)_?*_?*_checkpoint.pvtu $(FLMLMODEL1)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL1)_?*_checkpoint.ele $(FLMLMODEL1)_?*_checkpoint.face $(FLMLMODEL1)_?*_checkpoint.node
	rm -f $(FLMLMODEL1)_?*_checkpoint.vtu
	rm -f $(FLMLMODEL2)_?*.pvtu $(FLMLMODEL2)_?*.vtu
	rm -f $(FLMLMODEL2).detectors $(FLMLMODEL2).stat tmpf*
	rm -f $(FLMLMODEL2).log
	rm -f $(FLMLMODEL2)_?*_checkpoint.flml
	rm -f $(FLMLMODEL2)_?*_checkpoint_?*.ele $(FLMLMODEL2)_?*_checkpoint_?*.face $(FLMLMODEL2)_?*_checkpoint_?*.halo $(FLMLMODEL2)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL2)_?*_?*_checkpoint.pvtu $(FLMLMODEL2)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL2)_?*_checkpoint.ele $(FLMLMODEL2)_?*_checkpoint.face $(FLMLMODEL2)_?*_checkpoint.node
	rm -f $(FLMLMODEL2)_?*_checkpoint.vtu
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.face adapted_mesh.ele adapted_mesh.node
	rm -f adapted_mesh_?*.face adapted_mesh_?*.ele adapted_mesh_?*.node
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
