include envcheck.mk

FILENAME=MonaiValley_C_p1p1_nu0.01_kmkstab_drag0.002_butcircularoundisland0.2

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 2d:
	cd mesh; gmsh -bin MonaiValley_C.geo -2; ${BINPREFIX}fldecomp -n 4 -m gmsh MonaiValley_C

run: envcheck
	@echo **********Calling fluidity with verbose log output enabled:
	mpirun -n 4 ${BINPREFIX}fluidity -v2 -l MonaiValley_C_p1p1_nu0.01_kmkstab_drag0.002_butcircularoundisland0.2.flml

postprocess: envcheck
	@echo **********Calling the data extraction and plotting python script:
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH python ./plotgages.py --file=MonaiValley_C_p1p1_nu0.01_kmkstab_drag0.002_butcircularoundisland0.2.detectors

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *.pdf *.png *.stat *.vtu *checkpoint* *.log-* *.err-* matrixdump* *detectors *dat $(FILENAME)_*
	cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.halo
	@echo **********Cleaning the output from previous python runs:
	rm -rf *.pyc

