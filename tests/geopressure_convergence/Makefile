default: input

input:
	tar -xvjf src/meshes.tar.bz2
	../../bin/triangle2gmsh square-8
	../../bin/triangle2gmsh square-16
	../../bin/triangle2gmsh square-32
	../../bin/triangle2gmsh square-64
	cp geopressure_8.flml geopressure_16.flml 
	sed -i "s/"geopressure_8"/"geopressure_16"/g" geopressure_16.flml 
	sed -i "s/square-8/square-16/g" geopressure_16.flml 
	cp geopressure_8.flml geopressure_32.flml 
	sed -i "s/"geopressure_8"/"geopressure_32"/g" geopressure_32.flml 
	sed -i "s/square-8/square-32/g" geopressure_32.flml 
	cp geopressure_8.flml geopressure_64.flml 
	sed -i "s/"geopressure_8"/"geopressure_64"/g" geopressure_64.flml 
	sed -i "s/square-8/square-64/g" geopressure_64.flml


clean: clean-mesh clean-run-debug
	rm -f geopressure_16.flml
	rm -f geopressure_32.flml
	rm -f geopressure_64.flml
	rm -rf *flredecomp*
clean-mesh:
	rm -f square-?*.edge square-?*.ele square-?*.node
	rm -f square-?*_?*.edge square-?*_*.ele square-?*_?*.node square-?*_?*.halo
clean-run:
	rm -f geopressure_8_?*_?*.vtu rm -f geopressure_8_?*.pvtu
	rm -f geopressure_8.detectors geopressure_8.stat tmpf*
	rm -f geopressure_16_?*_?*.vtu rm -f geopressure_16_?*.pvtu
	rm -f geopressure_16.detectors geopressure_16.stat tmpf*
	rm -f geopressure_32_?*_?*.vtu rm -f geopressure_32_?*.pvtu
	rm -f geopressure_32.detectors geopressure_32.stat tmpf*
	rm -f geopressure_64_?*_?*.vtu rm -f geopressure_64_?*.pvtu
	rm -f geopressure_64.detectors geopressure_64.stat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info
	rm -f *.log
	rm -rf geopressure_8_?
	rm -rf geopressure_16_?
	rm -rf geopressure_32_?
	rm -rf geopressure_64_?
clean-run-debug: clean-run
