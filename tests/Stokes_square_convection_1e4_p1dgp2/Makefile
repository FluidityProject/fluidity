FLFLAGS = -v2 -l
FLMLMODEL = stokes-1e4-p1dgp2
FLUIDITY = dfluidity
MESH = square-structured-linear
MODEL = Stokes-square-convection-1e4-p1dgp2

default: input

input: src/$(MESH).msh clean
	cp src/$(MESH).msh .

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).msh
clean-run:
	rm -f $(FLMLMODEL).detectors $(FLMLMODEL)_?*.vtu $(FLMLMODEL)_?*_checkpoint.msh $(FLMLMODEL)_?*_checkpoint.flml $(FLMLMODEL)_?*_?*_checkpoint.vtu $(FLMLMODEL).stat fluidity.err-0 fluidity.log-0 matrixdump matrixdump.info
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu

run: clean-run input
	$(FLUIDITY) $(FLFLAGS) $(MODEL).flml
run-debug: clean-run-debug input
	$(FLUIDITY)-debug $(FLFLAGS) $(MODEL).flml
