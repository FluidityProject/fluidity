input:	clean
	./copy_script
	gmsh -2 -bin src/MMS_A.geo
	gmsh -2 -bin src/MMS_B.geo
	gmsh -2 -bin src/MMS_C.geo

run:
	fluidity MMS_A.flml
	fluidity MMS_B.flml
	fluidity MMS_C.flml

clean:
	rm -f *.vtu *.stat *.s *_A* *_B* *_C* *_D* src/*A.* src/*B.* src/*C.* src/*D.* *.log matrix* *~ *.pyc fluidity.* nohup.out *.node *.ele *.edge
