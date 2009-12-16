PROJECT = gyre

input:	$(PROJECT).ele

$(PROJECT).ele: $(PROJECT).msh
	../../scripts/gmsh2triangle --2d $<

$(PROJECT).msh: src/$(PROJECT).geo
	gmsh -2 -o $@ $<

clean:
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.edge *.node *.msh
