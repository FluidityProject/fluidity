FLMLMODEL = Eddying_TestCase

default: input

input: clean
	tar -xvjf src/InputCheckpoint.tar.bz2

clean: clean-checkpoint clean-run-debug
clean-checkpoint:
	rm -f Eddying_Detect_checkpoint_checkpoint_*
clean-run:
	rm -f $(FLMLMODEL)_?*.{p,}vtu
	rm -f $(FLMLMODEL).{detectors,stat} tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f $(FLMLMODEL)_?*_checkpoint.flml
	rm -f $(FLMLMODEL)_?*_checkpoint_?*.{msh,halo}
	rm -f $(FLMLMODEL)_?*_?*_checkpoint.{pvtu,vtu}
	rm -f $(FLMLMODEL)_?*_checkpoint.{msh,halo}
	rm -f $(FLMLMODEL)_?*_checkpoint.vtu
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu
