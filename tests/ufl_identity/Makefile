BACKEND ?= sequential

input: clean
	
	gmsh -2 -clscale 1.0 -o identity.msh src/cdisk.geo 
	../../bin/gmsh2triangle --2d identity.msh
	sed -e 's/{SIMULATION}/ufl_identity/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/1/g' src/identity.flml.template > identity.flml

clean:
	rm -f *.flml *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.msh
