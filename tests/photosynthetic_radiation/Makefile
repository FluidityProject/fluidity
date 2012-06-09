input: clean
	../../bin/gmsh2triangle src/column.msh
	../../bin/fldecomp -n 2 column

clean:
	rm -f *.face *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat
