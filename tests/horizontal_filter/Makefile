input: clean
	cd src; gmsh -2 square.geo ; cd ..;
	cp src/square.msh .

clean:
	rm -f *.ele *.edge *.node *.msh *.poly *.vtu *.s *.d.1 *.stat	fluidity.[le]* \
	matrixdump matrixdump.info
