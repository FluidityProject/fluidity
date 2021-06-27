input: default
default: clean
	gmsh -3 src/cylinder.geo
	cp src/cylinder.msh .

clean: 
	rm -f *.msh *.vtu *.log *.stat
