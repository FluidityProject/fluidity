input:
	cp src/*.msh src/*.halo .

inputbox:
	cp src/box_*.msh src/box_*.halo .
inputushape:
	cp src/ushape_*.msh src/ushape_*.halo .
inputucircle:
	cp src/circle.msh .

clean:clean-light
	rm -rf *.vtu *.pvtu *.stat

clean-light:
	rm -rf *.msh *.halo *.exo *.log* *.err* box_*_VolumeFraction ushape_*_VolumeFraction circle_*_VolumeFraction

test1:clean-light inputbox inputushape
	mpiexec -n 2 ../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m gmsh box ushape
	mv diagnostics.stat diag_test1.stat;

test2:clean-light inputbox inputushape
	mpiexec -n 2 ../../bin/nonmatching_domain_interpolation -v 3 -l -f 0 -m gmsh box ushape
	mv diagnostics.stat diag_test2.stat;

test3:clean-light inputbox inputucircle
	mpiexec -n 2 ../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m gmsh box circle
	mv diagnostics.stat diag_test3.stat;

test4:clean-light inputbox inputucircle
	mpiexec -n 2 ../../bin/nonmatching_domain_interpolation -v 3 -l -f 0 -m gmsh box circle
	mv diagnostics.stat diag_test4.stat;



test5:clean-light inputbox inputushape
	cp src/rendezvous-advfront-alg-test-fluidbox.exo src/ushape.exo .
	mpiexec -n 2 ../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m exodusii box ushape
	mv diagnostics.stat diag_test5.stat;

test6:clean-light inputbox inputushape
	mpiexec -n 2 ../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m gmsh ushape box
	mv diagnostics.stat diag_test6.stat;
