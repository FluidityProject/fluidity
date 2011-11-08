ASPECT=1

input: clean
	sed src/shelf3dbase.geo -e "s/^\(aspect *= *\)[0-9\.E\-][0-9\.E\-]*\(;\)/\1${ASPECT}\2/" > ./src/shelf3d.geo
	gmsh -3 ./src/shelf3d.geo
	../../bin/gmsh_mesh_transform "shelflength = 550, aspect = ${ASPECT}, shelfslopeheight = 800, minoceandepth = 100" 'x < (shelflength * aspect)' '(x, y, (z/(shelfslopeheight + minoceandepth + 200)) * ((x/(shelflength * aspect)) * shelfslopeheight + minoceandepth))' ./src/shelf3d.msh
	../../bin/gmsh2triangle ./src/shelf3d.msh

clean: 
	rm -f  ./src/*.{msh,bak} ./vtu/*.{,p}vtu *.ele *.edge *.node *.face *.poly *.{p,}vtu *.s  *.stat *.log-0 *.err-0 matrixdump*

