input:	clean
	gmsh -2 -bin src/MMS_A.geo
	gmsh -2 -bin src/MMS_B.geo
	gmsh -2 -bin src/MMS_C.geo

run:
	../../bin/fluidity MMS_A.flml
	../../bin/fluidity MMS_B.flml
	../../bin/fluidity MMS_C.flml

clean:
	rm -f *.vtu *.stat *.s *A_* *B_* *C_* *.log matrix* *~ mms_les_second_order_p2p1.py *.pyc fluidity.* nohup.out
