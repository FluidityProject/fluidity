input: clean
	cd src; gmsh -2 mesh2.geo

clean:
	rm -rf *.d.* *.s* *.vtu *.log *.edge *.node *.ele *.face *.vtu matrixdump* *.pyc src/*.msh fluidity.* *.detectors
	rm -rf *checkpoint*.flml
	rm -rf src/*.ele src/*.node src/*.edge src/*.msh src/*.face *.bak

