input: clean
	gmsh -2 ./src/incline_2-5.geo
	../../scripts/gmsh2triangle --2d ./src/incline_2-5.msh
	python ./src/box2incline_2D.py

clean: 
	rm -f ./src/*.msh *.ele *.edge *.node *.face *.poly *.vtu *.s  *.stat *.log* *.err-0 matrixdump* *checkpoint*
