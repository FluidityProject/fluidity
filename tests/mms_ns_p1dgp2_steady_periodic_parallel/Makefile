input:	clean
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo
	gmsh -2 src/MMS_D.geo
	cp src/MMS_A.msh .
	cp src/MMS_B.msh .
	cp src/MMS_C.msh .
	cp src/MMS_D.msh .
	periodise MMS_A.flml
	periodise MMS_B.flml
	periodise MMS_C.flml
	periodise MMS_D.flml

clean:
	rm -f *.vtu *.stat *.s *.log *.node *.ele *.edge *checkpoint*
