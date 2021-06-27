input:	clean
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo
	gmsh -2 src/MMS_D.geo
	# gmsh -2 src/MMS_E.geo
	cp src/MMS_A.msh .
	cp src/MMS_B.msh .
	cp src/MMS_C.msh .
	cp src/MMS_D.msh .
	# cp src/MMS_E.msh .

clean:
	rm -f *.vtu *.stat *.s *.log *.msh *checkpoint*
