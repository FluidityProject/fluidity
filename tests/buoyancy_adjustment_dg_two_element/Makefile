
input: clean
	gmsh -2 src/square.geo -o square.msh

clean:
	rm -f  *.ele *.edge *.node *.msh *.vtu *.stat matrixdump* *.node.bak fluidity.*

run: input
	fluidity -v3 -l mix.flml

