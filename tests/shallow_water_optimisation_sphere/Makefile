input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp sw_fsphere_A.swml sw_fsphere_B.swml; sed -i 's/_A/_B/' sw_fsphere_B.swml; spud-set sw_fsphere_B.swml /timestepping/timestep 2235.70667 
	cp sw_fsphere_A.swml sw_fsphere_C.swml; sed -i 's/_A/_C/' sw_fsphere_C.swml; spud-set sw_fsphere_C.swml /timestepping/timestep 1117.85334
	cp optimiser_A.oml optimiser_B.oml; sed -i 's/_A/_B/' optimiser_B.oml
	cp optimiser_A.oml optimiser_C.oml; sed -i 's/_A/_C/' optimiser_C.oml
clean:
	rm -f *.vtu *.stat *~*~ *.html *.pyc  *.pkl *.dat

clean_all:
	rm -f *.vtu *.stat *~*~ *.html *.pyc sw_fsphere_B.swml sw_fsphere_C.swml optimiser_B.oml optimiser_C.oml
	cd src; rm *ele *node *bound *~*~
