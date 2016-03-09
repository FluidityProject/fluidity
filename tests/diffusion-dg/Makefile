input: default
default: clean
	cp src/MMS_A.msh .
	cp src/MMS_B.msh .
	cp src/MMS_C.msh .
	cp src/MMS_D.msh .
	cp src/MMS_E.msh .
stretch: default
	./stretch.py -f MMS_A -x 100 -y 1
	./stretch.py -f MMS_D -x 100 -y 1
clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat *.log-0 *.err-0 matrixdump*

geos:
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo
	gmsh -2 src/MMS_D.geo
	gmsh -2 src/MMS_E.geo

setfiles:
	python set_input_files.py