CHECKPOINTEDMODEL = coarse-corner
MODEL = coarse-corner

default: input

input: clean-checkpoint
	cp -r ../$(MODEL)_*6_checkpoint* .
	./change_options

clean: clean-checkpoint clean-run-debug
clean-checkpoint:
	rm -f $(CHECKPOINTEDMODEL)_?*_checkpoint.msh $(CHECKPOINTEDMODEL)_?*_checkpoint.halo $(CHECKPOINTEDMODEL)_?*_checkpoint.flml $(CHECKPOINTEDMODEL)_?*_?*_checkpoint.vtu
clean-run:
	rm -f $(MODEL)_?*.pvtu $(MODEL)_?/$(MODEL)_?*.vtu
	rm -rf $(MODEL)_?
	rm -f $(MODEL).stat
	rm -f fluidity.err-?* fluidity.log-?* PI*
	rm -f $(MODEL)_?*_checkpoint_?*.{msh,halo} $(MODEL)_?*_checkpoint.flml $(MODEL)_?*_?*_checkpoint.pvtu $(MODEL)_?*_?*_checkpoint/$(MODEL)_?*_?*_checkpoint.vtu
	rm -rf $(MODEL)_?*_?_checkpoint*
	rm -f matrixdump matrixdump.info
clean-run-debug: clean-run
	rm -f bounding_box_?*_?*.{pvtu,vtu} final_metric_?*.{pvtu,vtu} gradation_metric_?*.{pvtu,vtu} interpolation_metric_final_?*.{pvtu,vtu} interpolation_metric_hessian_?*_?*.{pvtu,vtu} interpolation_metric_merge_?*_?*.{pvtu,vtu} interpolation_metric_metric_?*_?*.{pvtu,vtu} metric_input_?*_?*.{pvtu,vtu}
	rm -f adapted_mesh_?*.{pvtu,vtu} {pre,post}_adapt_halos_?*.{msh,halo} {pre,post}_adapt_{linear,pressure}_mesh_?*_?*.{msh,halo}
	rm -f gmon.out
