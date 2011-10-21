input: clean
	gmsh -3 sphere.msh src/sphere.geo
	gmsh -3 fluid_domain.msh src/fluid_domain.geo
	../../bin/gmsh2triangle sphere.msh
	../../bin/gmsh2triangle fluid_domain.msh

clean:
	rm -rf *.d.* *.s* *.vtu *.pvtu *.log *.halo *.node *.ele *.edge *.face matrixdump* *.pyc *.msh fluidity.* *.detectors *_data *checkpoint*
