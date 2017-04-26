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
	rm -f *.vtu *.stat *.s *_A* *_B* *_C* *_D* src/*A.* src/*B.* src/*C.* src/*D.* *.log matrix* *~ *.pyc fluidity.* nohup.out *.node *.ele *.edge
