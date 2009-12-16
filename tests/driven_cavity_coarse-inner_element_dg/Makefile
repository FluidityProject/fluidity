input: clean
	../../scripts/gmsh2triangle src/mesh.msh

clean:
	rm -rf *.d.* *.s* *.vtu *.log-0 *.err-0 mesh.node mesh.ele mesh.face 


# -*~ EXTRAS ~*-

run: input
	../../bin/dfluidity -l -v1 --xml driven_cavity_coarse-inner_element_dg.flml

test:
	cd ../../; time tools/testharness.py -f driven_cavity_coarse-inner_element_dg.xml --just-test; cd -

all: run
	cd ../../; time tools/testharness.py -f driven_cavity_coarse-inner_element_dg.xml --just-test; cd -

