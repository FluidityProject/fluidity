GMSH_GEO_FILE=src/spherical_segment.geo
NPROCS=3
input: clean
	gmsh -2 -bin $(GMSH_GEO_FILE) -o patch.msh
	../../bin/fldecomp -n $(NPROCS) -s -m gmsh patch

clean:
	rm -rf *.msh *.halo *.vtu *.stat fluidity.*
	rm -rf spherical_segment_* *.pvtu
	rm -rf flredecomp.* matrixdump*

