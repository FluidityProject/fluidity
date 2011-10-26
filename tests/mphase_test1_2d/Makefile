
BIN = `which dfluidity`

mesh: distclean
	gmsh -2 src/test.geo
	../../bin/gmsh2triangle --2d src/test.msh

p1p1:
	${BIN} -v3 -x test_p1p1.flml | tee test_p1p1.log

anisotropic_p1p1:
	${BIN} -v3 -x test_anisotropic_p1p1.flml | tee test_anisotropic_p1p1.log

p1dgp2:
	${BIN} -v3 -x test_p1dgp2.flml | tee test_p1dgp2.log

clean:
	rm -rf *.stat *.vtu *.log matrix* *checkpoint* *.convergence *.log

distclean: clean
	rm -rf *.node *.edge *.ele src/*.msh

