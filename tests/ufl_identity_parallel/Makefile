BACKEND ?= sequential

input: clean
	gmsh -2 -clscale 1.0 -o identity.msh src/cdisk.geo 
	../../bin/gmsh2triangle --2d identity.msh
	../../bin/fldecomp -n 3 identity
	sed -e 's/{SIMULATION}/ufl_identity_parallel/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/1/g' src/identity.flml.template > identity.flml

clean:
	rm -f *.halo *.flml *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh
	rm -rf ufl_identity_parallel_[01] *.pvtu 
