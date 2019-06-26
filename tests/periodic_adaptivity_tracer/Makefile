input:
	gmsh -format msh2 -2 src/mesh.geo -o mesh.msh

clean:
	rm -f *.vtu *.log-? *.err-? *.stat *.msh *.edge *.ele *.node
	rm -f periodic_adaptivity_tracer_periodised.flml \
	matrixdump matrixdump.info
