default: input

input: src/target.msh
	../../bin/gmsh2triangle --2d src/target.msh

clean: clean-mesh clean-run
clean-mesh:
	rm -f target.edge target.ele target.node
clean-run:
	rm -f matrixdump matrixdump.info
	rm -f output.vtu
	rm -f tmpf*
	rm -f *.log
