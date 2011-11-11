PROJECT = gyre

input:	$(PROJECT).ele

$(PROJECT).ele: $(PROJECT).msh
	../../bin/gmsh2triangle $<

$(PROJECT).msh: src/$(PROJECT).geo
	gmsh -3 -o $@ $<

clean:
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.edge *.node *.msh
