input:
	cp src/*.msh .

inputbox:
	cp src/box_*.msh src/box_*.halo .
inputushape:
	cp src/ushape.msh .
inputucircle:
	cp src/circle.msh .

clean:
	rm -rf *.vtu *.pvtu *.stat matrixdump* *.msh *.exo *.log* *.err*

clean-light:
	rm -rf *.msh *.exo *.log* *.err*

test1:clean-light inputbox inputushape
	../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m gmsh box ushape
	mv diagnostics.stat diag_test1.stat;

test2:clean-light inputbox inputushape
	../../bin/nonmatching_domain_interpolation -v 3 -l -f 0 -m gmsh box ushape
	mv diagnostics.stat diag_test2.stat;

test3:clean-light inputbox inputucircle
	../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m gmsh box circle
	mv diagnostics.stat diag_test3.stat;

test4:clean-light inputbox inputucircle
	../../bin/nonmatching_domain_interpolation -v 3 -l -f 0 -m gmsh box circle
	mv diagnostics.stat diag_test4.stat;



test5:clean-light inputbox inputushape
	cp src/rendezvous-advfront-alg-test-fluidbox.exo src/ushape.exo .
	../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m exodusii box ushape
	mv diagnostics.stat diag_test5.stat;

test6:clean-light inputbox inputushape
	../../bin/nonmatching_domain_interpolation -v 3 -l -f 1 -m gmsh ushape box
	mv diagnostics.stat diag_test6.stat;
