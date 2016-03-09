FLFLAGS = -v2 -l
FLMLMODEL = square-convection
FLUIDITY = dfluidity
MESH = square-structured-linear
MODEL = square-convection

default: input

input: clean
	cp src/* .

clean: clean-mesh clean-run-debug 
clean-mesh:
	rm -f $(MESH)_?*.edge $(MESH)_?*.ele $(MESH)_?*.node $(MESH)_?*.halo $(MESH)_?*.msh
clean-run:
	rm -f *.dat *.halo *.pvtu *.vtu *.edge *.ele *.node $(FLMLMODEL)_?*_checkpoint.flml *.detectors *.stat *.err-? *.log-? matrixdump matrixdump.info
	rm -rf $(MODEL)*_checkpoint $(MODEL)_?
	rm -rf *flredecomp*
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu

run: clean-run input
	$(FLUIDITY) $(FLFLAGS) $(MODEL).flml
run-debug: clean-run-debug input
	$(FLUIDITY)-debug $(FLFLAGS) $(MODEL).flml
