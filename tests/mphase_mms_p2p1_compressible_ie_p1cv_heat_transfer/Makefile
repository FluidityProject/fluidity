input:	clean
	./copy_script
	gmsh -2 -bin src/MMS_A.geo
	gmsh -2 -bin src/MMS_B.geo
	gmsh -2 -bin src/MMS_C.geo

run:
	@echo **********Running simulations
	mpiexec fluidity MMS_A.flml
	mpiexec fluidity MMS_B.flml
	mpiexec fluidity MMS_C.flml

clean:
	rm -f *.vtu *.stat *.s *.log *_A* *_B* *_C* matrix* *~ mphase_mms_p2p1_compressible_ie_p1cv_heat_transfer.py *.pyc fluidity.*
