input: default
default: clean
	cd src
	gmsh -2 square.geo
	cd ..

clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat
	rm -f *.err-? *.log-? \
	matrixdump matrixdump.info
