PROJECT = checkmesh

input: clean

	gmsh -3 src/good_mesh.geo -o good_mesh.msh
	cp good_mesh.msh bad_mesh.msh
	../../bin/gmsh_mesh_transform 'True' '(x,x*y,x*z)' bad_mesh.msh 

clean:

	rm -f *.msh *.geo *.bak *.log
