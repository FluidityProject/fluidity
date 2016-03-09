input: clean
	
checkpoint_input: clean
	gmsh -3 -o gls-StationPapa-3D.msh src/gls-StationPapa-3D.geo

clean: 
	rm -f gls-StationPapa-3D.* gls-StationPapa-param_*.vtu gls-StationPapa-orig_*.vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *log *err
