input: clean
	gmsh -format msh2 -2 -bin src/square.geo -o square.msh

clean:
	rm -rf output* square.msh pvtu2vtu_flredecomp*
	rm -rf flredecomp.* fluidity.*
