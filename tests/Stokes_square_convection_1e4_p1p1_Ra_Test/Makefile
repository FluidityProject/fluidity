FLFLAGS = -v2 -l
FLMLMODEL = stokes-1e4-p1p1
FLUIDITY = fluidity
MESH = square-structured-linear
MODEL = Stokes-square-convection-1e4-p1p1-Ra-Test

default: input

input: src/$(MESH).edge src/$(MESH).ele src/$(MESH).node clean
	cp src/$(MESH).edge src/$(MESH).ele src/$(MESH).node .
	triangle2gmsh $(MESH)

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).edge $(MESH).ele $(MESH).node $(MESH).msh
clean-run:
	rm -f $(FLMLMODEL)_?*.vtu $(FLMLMODEL)_?*_checkpoint.msh $(FLMLMODEL)_?*_checkpoint.flml $(FLMLMODEL)_?*_?*_checkpoint.vtu $(FLMLMODEL).stat fluidity.err-0 fluidity.log-0 matrixdump matrixdump.info *.vtu *.stat *.detectors
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu

run: clean-run input
	$(FLUIDITY) $(FLFLAGS) $(MODEL).flml
run-debug: clean-run-debug input
	$(FLUIDITY)-debug $(FLFLAGS) $(MODEL).flml
