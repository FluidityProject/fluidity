FC = $(shell cat ../../Makefile | grep ^FC | awk '{print $$3}')

input: clean
	$(FC) src/set_phi.f -o src/set_phi
	gmsh -3 src/mmat-gravity-col.geo
	../../scripts/gmsh2triangle src/mmat-gravity-col.msh

clean:
	rm -rf *.d.* *.stat *.vtu *.convergence mmat-gravity-col.node mmat-gravity-col.face mmat-gravity-col.ele src/set_phi *checkpoint*
