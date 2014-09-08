input: clean
	gmsh -2 -clmax 0.2 box2d.geo -o box2dA.msh
	gmsh -2 -clmax 0.1 box2d.geo -o box2dB.msh
	gmsh -2 -clmax 0.05 box2d.geo -o box2dC.msh
	gmsh -2 -clmax 0.025 box2d.geo -o box2dD.msh
	gmsh -2 -clmax 0.0125 box2d.geo -o box2dE.msh	

clean:
	rm -f fluidity.*
	rm -f  *.stat *.vtu
	rm -f  matrixdump*

