input: clean
	gmsh -2 ./src/lock_exchange_2d.geo
	cp src/lock_exchange_2d.msh .

clean: 
	rm -f ./src/*.msh *.msh *.poly *.vtu *.s  *.stat *.log* *.err-0 matrixdump* *checkpoint*
