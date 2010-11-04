nprocs=64
preprocess: 
#	gmsh -2 med.geo
	../../scripts/gmsh2triangle med.msh
	sed < med.node 's/[0-9]*/2/2' > med2.node
	mv med2.node med.node
ifeq ($(nprocs),1)
	@echo **********Serial run: not decomposing mesh
else
	@if [ ! -e ../../bin/fldecomp ] ; then echo "You must build Fluidity parallel tools, using the command 'make fltools' in the Fluidity directory, prior to running this test in parallel."; false; fi
	@echo **********Decomposing the mesh into $(nprocs) parts for parallel run:
	../../bin/fldecomp -n $(nprocs) -s med 
endif

run:
ifeq ($(nprocs),1)
	@echo **********Calling fluidity in serial with verbose log output enabled:
	../../bin/fluidity -v2 -l med.flml
else
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpirun -np $(nprocs)  ../../bin/fluidity -v2 -l med.flml
endif

postprocess:
	@echo **********Reprojecting final output to long-lat and extracting RMS differences:
	python ./reproject-to-loglat.py
	python ./Med-tides-probe.py


clean:
	rm -rf  *gmt* *.ele *.node *.edge Med-GEBCO-5m-gauges-fes2004-O1-102 *med_* *.halo *.pvtu *.vtu *.stat *.log-* *.err-* *check* *Check* *matrixdump*

input: clean preprocess
