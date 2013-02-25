input:	clean
	./copy_script
	gmsh -2 -bin src/MMS_A.geo
	gmsh -2 -bin src/MMS_B.geo
	gmsh -2 -bin src/MMS_C.geo

run:
	../../bin/fluidity MMS_A.flml
	../../bin/fluidity MMS_B.flml
	../../bin/fluidity MMS_C.flml

clean:
	rm -f *.vtu *.stat *.s *.log *_A* *_B* *_C* src/*A.* src/*B.* src/*C.* matrix* *~ mms_rans_p2p1_stress_form.py *.pyc fluidity.*
