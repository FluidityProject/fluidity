input:	clean
	./copy_script
	gmsh -2 -bin src/MMS_A.geo
	gmsh -2 -bin src/MMS_B.geo
	gmsh -2 -bin src/MMS_C.geo

clean:
	rm -f *.vtu *.stat *.s *.log *_A* *_B* *_C* *_D* matrix* *~ *.pyc fluidity.*

run:
	../../bin/fluidity MMS_A.flml
	../../bin/fluidity MMS_B.flml
	../../bin/fluidity MMS_C.flml
