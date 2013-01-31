DEGREE ?= 1
BACKEND ?= sequential

.PHONY: cdisk clean flml unitsquare

clean:
	rm -f *.halo *.flml *.ele *.edge *.node *.poly *.pvtu *.vtu *.s *.d.1 *.stat *.msh *.geo *.pyc
	rm -rf MMS_[ABCD]_[0-9]

cdisk:
	gmsh -2 -clscale 1.0  -o MMS_D.msh ufl_common/cdisk.geo
	gmsh -2 -clscale 2.0  -o MMS_C.msh ufl_common/cdisk.geo
	gmsh -2 -clscale 4.0  -o MMS_B.msh ufl_common/cdisk.geo
	gmsh -2 -clscale 8.0 -o MMS_A.msh ufl_common/cdisk.geo
	../../bin/gmsh2triangle --2d MMS_A.msh
	../../bin/gmsh2triangle --2d MMS_B.msh
	../../bin/gmsh2triangle --2d MMS_C.msh
	../../bin/gmsh2triangle --2d MMS_D.msh

parallel-cdisk: cdisk
	../../bin/fldecomp -n 3 MMS_A
	../../bin/fldecomp -n 3 MMS_B
	../../bin/fldecomp -n 3 MMS_C
	../../bin/fldecomp -n 3 MMS_D

flml:
	sed -e 's/{MESH}/MMS_A/g' -e 's/{SIMULATION}/MMS_A/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/$(DEGREE)/g' ufl_common/$(FLML) > MMS_A.flml
	sed -e 's/{MESH}/MMS_B/g' -e 's/{SIMULATION}/MMS_B/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/$(DEGREE)/g' ufl_common/$(FLML) > MMS_B.flml
	sed -e 's/{MESH}/MMS_C/g' -e 's/{SIMULATION}/MMS_C/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/$(DEGREE)/g' ufl_common/$(FLML) > MMS_C.flml
	sed -e 's/{MESH}/MMS_D/g' -e 's/{SIMULATION}/MMS_D/g' -e 's/{BACKEND}/$(BACKEND)/g' -e 's/{DEGREE}/$(DEGREE)/g' ufl_common/$(FLML) > MMS_D.flml

unitsquare:
	ufl_common/generate_mesh MMS_A 30
	ufl_common/generate_mesh MMS_B 60
	ufl_common/generate_mesh MMS_C 120
	ufl_common/generate_mesh MMS_D 240
