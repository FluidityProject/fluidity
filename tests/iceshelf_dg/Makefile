input: clean
	gmsh -2 ./src/shelf2d.geo
	../../bin/gmsh_mesh_transform 'shelflength = 550, shelfslopeheight = 800, minoceandepth = 100' 'x < shelflength' '(x, (y/(shelfslopeheight + minoceandepth)) * ((x/shelflength) * shelfslopeheight + minoceandepth) ,z)' ./src/shelf2d.msh
	../../bin/gmsh2triangle --2d ./src/shelf2d.msh

clean: 
	rm -f  ./src/*.{msh,bak} ./vtu/*.{,p}vtu *.ele *.edge *.node *.face *.poly *.{p,}vtu *.s  *.stat *.log-0 *.err-0 matrixdump*



# #ROOT = /d/src/fluidity
# ROOT = ${HOME}/cur/fluidity_iceshelf
# 
# default: input
# 
# input: clean
# 	gmsh -2 ./src/shelf2d.geo
# 	${ROOT}/scripts/gmsh_mesh_transform 'shelflength = 550, shelfslopeheight = 800, minoceandepth = 100' 'x < shelflength' '(x, (y/(shelfslopeheight + minoceandepth)) * ((x/shelflength) * shelfslopeheight + minoceandepth) ,z)' ./src/shelf2d.msh
# 	${ROOT}/scripts/gmsh2triangle --2d ./src/shelf2d.msh
# 
# 	#gmsh -3 -optimize -o icecube_no_flux.msh geos/ice_simple3D_3walls3.geo 
# 	#${ROOT}/scripts/gmsh2triangle icecube_no_flux.msh
# 
# clean: 
# 	rm -f  ./src/*.{msh,bak} ./vtu/*.{,p}vtu *.ele *.edge *.node *.face *.poly *.{p,}vtu *.s  *.stat *.log-0 *.err-0 matrixdump*
# 
# run: input
# 	mkdir -p ./vtu
# 	${ROOT}/bin/fluidity -v3 -l shelf2d.flml 
# 
# 	#${ROOT}/bin/fluidity -v2 -l test_caseCiceGravity_in.flml
# 
# diamond:
# 	diamond -s ${ROOT}/schemas/fluidity_options.rng shelf2d.flml
# 
# #/d/srcbzr/myfluidity_icer/bin/fluidity -v2 -l test_caseCiceGravity.flml
# 
# #	../../bin/gmsh_mesh_transform 'shelflength = 550, shelfslopeheight = 800, minoceandepth = 100' 'x < shelflength' '(x, y, (z/(shelfslopeheight + minoceandepth)) * ((x/shelflength) * shelfslopeheight + minoceandepth))' ./src/shelf2d.msh
#
