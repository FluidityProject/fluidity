input: clean
	gmsh -3 src/region_id.geo
	cp src/region_id.msh .

clean:
	rm -rf *.d.* *.s* *.vtu region_id.node region_id.ele region_id.face src/region_id.msh fluidity.[le]* \
	matrixdump matrixdump.info
