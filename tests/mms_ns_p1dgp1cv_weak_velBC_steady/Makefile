input:	clean
	../../bin/gmsh2triangle --2d src/MMS_A.msh
	../../bin/gmsh2triangle --2d src/MMS_B.msh
	../../bin/gmsh2triangle --2d src/MMS_C.msh

clean:
	rm -f *.vtu *.stat *.log *.node *.ele *.edge *checkpoint*
	rm -f *.msh
	rm -f matrixdump*
	rm -f fluidity.*
