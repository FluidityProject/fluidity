FLMLMODEL1 = geopressure_8
FLMLMODEL2 = geopressure_16
FLMLMODEL3 = geopressure_32
FLMLMODEL4 = geopressure_64

default: input

input:
	tar -xvjf src/meshes.tar.bz2
	../../bin/triangle2gmsh square-8
	../../bin/triangle2gmsh square-16
	../../bin/triangle2gmsh square-32
	../../bin/triangle2gmsh square-64
	cp $(FLMLMODEL1).flml $(FLMLMODEL2).flml 
	sed -i "s/"$(FLMLMODEL1)"/"$(FLMLMODEL2)"/g" $(FLMLMODEL2).flml 
	sed -i "s/square-8/square-16/g" $(FLMLMODEL2).flml 
	cp $(FLMLMODEL1).flml $(FLMLMODEL3).flml 
	sed -i "s/"$(FLMLMODEL1)"/"$(FLMLMODEL3)"/g" $(FLMLMODEL3).flml 
	sed -i "s/square-8/square-32/g" $(FLMLMODEL3).flml 
	cp $(FLMLMODEL1).flml $(FLMLMODEL4).flml 
	sed -i "s/"$(FLMLMODEL1)"/"$(FLMLMODEL4)"/g" $(FLMLMODEL4).flml 
	sed -i "s/square-8/square-64/g" $(FLMLMODEL4).flml 

clean: clean-mesh clean-run-debug
	rm -f $(FLMLMODEL2).flml
	rm -f $(FLMLMODEL3).flml
	rm -f $(FLMLMODEL4).flml
clean-mesh:
	rm -f square-?*.edge square-?*.ele square-?*.node
clean-run:
	rm -f $(FLMLMODEL1)_?*.vtu
	rm -f $(FLMLMODEL1).detectors $(FLMLMODEL1).stat tmpf*
	rm -f $(FLMLMODEL2)_?*.vtu
	rm -f $(FLMLMODEL2).detectors $(FLMLMODEL2).stat tmpf*
	rm -f $(FLMLMODEL3)_?*.vtu
	rm -f $(FLMLMODEL3).detectors $(FLMLMODEL3).stat tmpf*
	rm -f $(FLMLMODEL4)_?*.vtu
	rm -f $(FLMLMODEL4).detectors $(FLMLMODEL4).stat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f *.log
clean-run-debug: clean-run
