include envcheck.mk

preprocess: envcheck
	@echo **********Converting the gmsh mesh to triangle format:

run: envcheck
	@echo **********Calling fluidity with verbose log output enabled:
	${BINPREFIX}fluidity -v2 -l flow_past_sphere_Re1.flml

postprocess:
	@echo **********Calling the data extraction and plotting script
	python ./plot_data.py

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *.vtu *.stat *.log-0 *.err-0 matrixdump* 
	rm -rf first_timestep_adapted_mesh*
	rm -f Sphere_drag.pdf
