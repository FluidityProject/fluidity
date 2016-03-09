input:	clean
	gmsh -2 -o MMS_A.msh src/MMS_A.geo
	gmsh -2 -o MMS_B.msh src/MMS_B.geo
	gmsh -2 -o MMS_C.msh src/MMS_C.geo
	gmsh -2 -o MMS_D.msh src/MMS_D.geo

clean:
	rm -f *.vtu *.stat *.s *.log *.node *.ele *.edge *.msh *checkpoint*
	rm -f *.msh
	rm -f *.err-? *.log-? \
	matrixdump matrixdump.info
