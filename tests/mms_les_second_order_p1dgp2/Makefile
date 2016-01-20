input:	clean
	./copy_script
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo
#	gmsh -2 src/MMS_D.geo

run:  input
	../../bin/fluidity -v2 MMS_A.flml > MMS_A.log &
	../../bin/fluidity -v2 MMS_B.flml > MMS_B.log &
	../../bin/fluidity -v2 MMS_C.flml > MMS_C.log &
#	../../bin/fluidity -v2 MMS_D.flml > MMS_D.log &

clean:
	rm -f *.vtu *.stat *.s *_A* *_B* *_C* *_D* src/*A.* src/*B.* src/*C.* src/*D.* *.log matrix* *~ mms_les_second_order_p1dgp2.py *.pyc fluidity.* nohup.out *.node *.ele *.edge
