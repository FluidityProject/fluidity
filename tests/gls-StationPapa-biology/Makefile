input: clean
	gmsh -3 -o gls-StationPapa-3D.msh src/gls-StationPapa-3D.geo
	../../scripts/gmsh2triangle gls-StationPapa-3D.msh

clean: 
	rm -f  *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *log *err
