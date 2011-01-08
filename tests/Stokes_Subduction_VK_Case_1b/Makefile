input: clean
	cp src/Subduction_Mesh.msh ./
	../../scripts/gmsh2triangle --2d -i Subduction_Mesh.msh

clean:
	rm -f *.ele *.edge *.node *.vtu *.stat *.msh *.detectors *.dat fluidity.* *_state matrixdump*
