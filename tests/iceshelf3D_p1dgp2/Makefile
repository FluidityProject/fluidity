
input: clean

	gmsh -3 -optimize -o ./src/shelf3d.msh ./src/sym_shelf3d.geo   

	../../bin/gmsh_mesh_transform 'shelflength = 550, shelfslopeheight = 800, minoceandepth = 100' 'x < shelflength' '(x, y,(z/(shelfslopeheight + minoceandepth)) * ((x/shelflength) * shelfslopeheight + minoceandepth))' ./src/shelf3d.msh
	../../bin/gmsh2triangle ./src/shelf3d.msh

clean: 
	rm -f  ./src/*.{msh,bak} ./vtu/*.{,p}vtu *.ele *.edge *.node *.face *.poly *.{p,}vtu *.s  *.stat *.log-0 *.err-0 matrixdump* *.halo *.log-* *.err-* *.vtu
