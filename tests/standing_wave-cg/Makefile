PROJECT = standing_wave

input: clean
	gmsh -3 -bin src/standing_wave.geo -o standing_wave.msh

clean:
	rm -f $(PROJECT) $(PROJECT)0 $(PROJECT).fld
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu *.detectors.h5part
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node *.msh \
	matrixdump matrixdump.info
