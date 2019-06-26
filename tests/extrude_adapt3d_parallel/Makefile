input: clean
	gmsh -format msh2 -2 src/square.geo -o square.msh

clean:
	rm -fr *.msh *.halo *.vtu *.pvtu
	rm -rf fluidity.* *.stat
	rm -rf matrixdump matrixdump.info
	rm -rf *flredecomp* *checkpoint*
	rm -rf box_*
