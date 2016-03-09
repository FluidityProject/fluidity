FLMLMODEL = annulus
MESH = annulus-structured-tanh-23x64

default: input

input: src/$(MESH).tar.bz2
	tar -xvjf src/$(MESH).tar.bz2
	../../bin/triangle2gmsh $(MESH)

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).ele $(MESH).face $(MESH).node
	rm -f $(MESH)_?*.ele $(MESH)_?*.face $(MESH)_?*.halo $(MESH)_?*.node decomp.log decomp_err.log
clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*.vtu
	rm -f $(FLMLMODEL).detectors $(FLMLMODEL).stat $(FLMLMODEL).steady_state $(FLMLMODEL).steady_state.dat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f $(FLMLMODEL)_?*_checkpoint.flml
	rm -f $(FLMLMODEL)_?*_checkpoint_?*.ele $(FLMLMODEL)_?*_checkpoint_?*.face $(FLMLMODEL)_?*_checkpoint_?*.halo $(FLMLMODEL)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL)_?*_?*_checkpoint.pvtu $(FLMLMODEL)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL)_?*_checkpoint.ele $(FLMLMODEL)_?*_checkpoint.face $(FLMLMODEL)_?*_checkpoint.node
	rm -f $(FLMLMODEL)_?*_checkpoint.vtu
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.face adapted_mesh.ele adapted_mesh.node
	rm -f adapted_mesh_?*.face adapted_mesh_?*.ele adapted_mesh_?*.node
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
