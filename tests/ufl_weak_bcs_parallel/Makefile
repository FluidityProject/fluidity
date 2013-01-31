BACKEND ?= sequential

input: clean
	gmsh -2 -clscale 1.0 -o weak_bcs.msh src/cdisk.geo 
	../../bin/gmsh2triangle --2d weak_bcs.msh
	../../bin/fldecomp -n 3 weak_bcs
	sed -e 's/{SIMULATION}/ufl_weak_bcs/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/1/g' src/weak_bcs.flml.template > weak_bcs.flml

clean:
	rm -f *.flml *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh
	rm -rf *.halo *.pvtu ufl_weak_bcs_[012]
