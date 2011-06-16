input: meshes copy_tests clean

meshes:
	cd src; make

copy_tests:
	cp sw_fsphere_A.swml sw_fsphere_B.swml; sed -i 's/_A/_B/' sw_fsphere_B.swml; spud-set sw_fsphere_B.swml /timestepping/timestep 0.125
	cp sw_fsphere_A.swml sw_fsphere_C.swml; sed -i 's/_A/_C/' sw_fsphere_C.swml; spud-set sw_fsphere_C.swml /timestepping/timestep 0.0625

clean:
	rm -f *.vtu *.stat *~*~ *.html *.pyc 

clean_all:
	rm -f *.vtu *.stat *~*~ *.html *.pyc sw_fsphere_B.swml sw_fsphere_C.swml 
	cd src; rm *ele *node *bound *~*~
