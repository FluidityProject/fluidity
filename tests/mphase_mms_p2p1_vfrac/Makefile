input:	clean
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo
	./copy_script

run:
	@echo **********Running simulations
	fluidity MMS_A.flml
	fluidity MMS_B.flml
	fluidity MMS_C.flml

clean:
	rm -f *.vtu *.stat *.s *.log *_A* *_B* *_C* matrix* *~ mphase_mms_p2p1_vfrac.py *.pyc fluidity.*
