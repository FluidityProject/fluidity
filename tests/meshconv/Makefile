PROJECT = checkmesh

input: clean

	mkdir -p exodus2triangle
	cp inputs/block.e exodus2triangle/
	mkdir -p exodus2gmsh
	cp inputs/block.e exodus2gmsh/
	mkdir -p gmsh2triangle
	cp inputs/block.msh gmsh2triangle/

input-clean:
	rm exodus2triangle/block.e
	rm exodus2gmsh/block.e
	rm gmsh2triangle/block.msh

clean:

	rm -f exodus2triangle/* exodus2gmsh/* gmsh2triangle/* 
	rm -fd exodus2triangle exodus2gmsh gmsh2triangle
