FLMLMODEL = prescribed_adaptivity

default: input

input:
	tar -xvjf src/meshes.tar.bz2
	for i in `seq 0 9`; do ../../bin/triangle2gmsh $$i ; done

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f *.ele *.edge *.node *.msh
clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*.vtu
	rm -f $(FLMLMODEL).detectors $(FLMLMODEL).stat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f $(FLMLMODEL)_?*_checkpoint.flml
	rm -f $(FLMLMODEL)_?*_checkpoint_?*.ele $(FLMLMODEL)_?*_checkpoint_?*.edge $(FLMLMODEL)_?*_checkpoint_?*.halo $(FLMLMODEL)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL)_?*_?*_checkpoint.pvtu $(FLMLMODEL)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL)_?*_checkpoint.ele $(FLMLMODEL)_?*_checkpoint.edge $(FLMLMODEL)_?*_checkpoint.node
	rm -f $(FLMLMODEL)_?*_checkpoint.vtu
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f adapted_mesh.edge adapted_mesh.ele adapted_mesh.node
	rm -f adapted_mesh_?*.edge adapted_mesh_?*.ele adapted_mesh_?*.node
	rm -f adapted_state_?*.pvtu adapted_state_?*.vtu
