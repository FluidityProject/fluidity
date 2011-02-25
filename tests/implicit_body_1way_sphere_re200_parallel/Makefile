input: clean
	gmsh -3 sphere.msh src/sphere.geo
	gmsh -3 fluid_domain.msh src/fluid_domain.geo
	../../scripts/gmsh2triangle sphere.msh
	../../scripts/gmsh2triangle fluid_domain.msh
	../../bin/fldecomp -n 2 -f fluid_domain

clean:
	rm -rf *.d.* *.s* *.vtu *.pvtu *.log *.halo *.node *.ele *.edge *.face matrixdump* *.pyc *.msh fluidity.* *.detectors *_data *checkpoint*
