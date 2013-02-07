BACKEND ?= sequential

input: clean
	
	gmsh -2 -clscale 1.0 -o laplace.msh src/cdisk.geo 
	../../bin/gmsh2triangle --2d laplace.msh
	sed -e 's/{SIMULATION}/ufl_laplace/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/1/g' src/laplace.flml.template > laplace.flml

clean:
	rm -f *.flml *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh
