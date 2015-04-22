input: clean
	gmsh -2 src/sloshing_tank_mesh1.geo
	cp src/sloshing_tank_mesh1.msh .

clean:
	rm -rf *.stat* *.vtu sloshing_tank_mesh1.node sloshing_tank_mesh1.ele sloshing_tank_mesh1.edge src/sloshing_tank_mesh1.msh *.detectors matrixdump* fluidity.log-0 fluidity.err-0
