input:
	gmsh -2 src/MMS_B.geo -o MMS_B.msh
	../../scripts/gmsh2triangle --2d MMS_B.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.msh
