preprocess: 
	@echo **********Generating the mesh using gmsh in 2d:
	gmsh -2 channel_2d.geo
	@echo **********Converting the gmsh mesh to triangle format:
	../../scripts/gmsh2triangle -2 channel_2d.msh 

run:
	@echo **********Calling fluidity with verbose log output enabled:
	../../bin/fluidity -v2 -l lock_exchange.flml

postprocess:
	@echo **********Calling the data extraction and plotting python script:
	python ./plot_data.py

input: preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.msh *.node *.ele *.edge lock_exchange_*.vtu *checkpoint* *convergence* *.log-* *.err-* matrixdump*
	@echo **********Cleaning the output from previous python runs:
	rm -rf diagnostics le_tools.pyc
