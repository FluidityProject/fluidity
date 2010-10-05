PROJECT = gyre

input:
	gmsh -2 -bin src/$(PROJECT).geo -o $(PROJECT).msh

clean:
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.edge *.node *.msh
