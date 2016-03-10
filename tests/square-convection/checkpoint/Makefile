CHECKPOINTEDMODEL = square-convection
FLMLMODEL = square-convection-picked-up

default: input

input: clean-checkpoint
	cp ../square-convection_*50_checkpoint* .
	./change_options

clean: clean-checkpoint clean-run-debug
clean-checkpoint:
	rm -f $(CHECKPOINTEDMODEL)_?*_checkpoint.edge $(CHECKPOINTEDMODEL)_?*_checkpoint.ele $(CHECKPOINTEDMODEL)_?*_checkpoint.node $(CHECKPOINTEDMODEL)_?*_checkpoint.flml $(CHECKPOINTEDMODEL)_?*_?*_checkpoint.vtu
clean-run:
	rm -f $(FLMLMODEL).detectors.dat $(FLMLMODEL)_?*.vtu $(FLMLMODEL)_?*_checkpoint.edge $(FLMLMODEL)_?*_checkpoint.ele $(FLMLMODEL)_?*_checkpoint.node $(FLMLMODEL)_?*_checkpoint.flml $(FLMLMODEL)_?*_?*_checkpoint.vtu $(FLMLMODEL).stat $(FLMLMODEL).detectors fluidity.err-0 fluidity.log-0 matrixdump matrixdump.info
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
	rm -f pre_adapt_linear_mesh_?*.ele pre_adapt_linear_mesh_?*.face pre_adapt_linear_mesh_?*.node
	rm -f post_adapt_linear_mesh_?*.ele post_adapt_linear_mesh_?*.face post_adapt_linear_mesh_?*.node
