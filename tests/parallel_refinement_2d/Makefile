input:
	gmsh -2 src/square.geo -o square.msh

clean:
	rm -rf fluidity.* square.msh refinep* refine_* refine.stat *.pvtu
	rm -rf matrixdump*
	rm -rf adapted_state_[0-9] gradation_metric_[0-9] bounding_box_[0-9] metric_input_1_[0-9]
