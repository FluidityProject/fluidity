NAME = porous_block_tri

MESH = ${PWD}/src/${NAME}.msh
FLML = ${PWD}/${NAME}.flml
FLOPTS = -v3
MESHOPTS =

mesh: clean
	../../scripts/gmsh2triangle ${MESHOPTS} ${MESH}

run: mesh
	../../bin/dfluidity ${FLOPTS} -x ${FLML} | tee ${NAME}.log

clean:
	rm -f  *.ele *.edge *.face *.node *.poly *.vtu *.s *.d.1 *.stat *.log

