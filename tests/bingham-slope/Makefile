input: clean
	gmsh -3 src/bingham-slope.geo
	../../scripts/gmsh2triangle src/bingham-slope.msh
	./box2incline.py

clean:
	rm -rf *.d.* *.s* *.vtu bingham-slope bingham-slope0 bingham-slope.fld bingham-slope.node bingham-slope.ele bingham-slope.face
