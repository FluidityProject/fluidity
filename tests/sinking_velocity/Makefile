input: clean
	gmsh -3 -bin -algo front3d -bin src/column.geo -o column.msh

clean:
	rm -f  *.ele *.edge *.face *.node *.poly *.vtu *.s *.d.1 *.stat *.msh \
	matrixdump matrixdump.info
