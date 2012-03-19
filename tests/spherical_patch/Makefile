GMSH_GEO_FILE=src/spherical_segment.geo
NPROCS=3
input:
	gmsh -2 -bin $(GMSH_GEO_FILE) -o patch.msh
	../../bin/fldecomp -n $(NPROCS) -s -m gmsh patch

clean:
	rm -rf *.msh *.halo

