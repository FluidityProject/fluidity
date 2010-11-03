preprocess:
	@echo **********Converting input meshes to triangle format for Re=1,10,100,1000
	$(foreach Re,1 10 100 1000, \
		cd sphere-3D-drag-Re$(Re); \
		../../../scripts/gmsh2triangle ../src/unit_sphere.msh; \
		cd ..; \
	)
	@echo **********Attempting to decompose the mesh into 8 partitions for the Re=1000 case
	cd sphere-3D-drag-Re1000; \
	../../../bin/fldecomp -n 8 -f unit_sphere; \
	cd ..; \

input: preprocess

postprocess:
	@echo **********No postprocessing needed


clean:
	$(foreach Re,1 10 100 1000, \
		cd sphere-3D-drag-Re$(Re); \
		rm -f *.ele *.edge *.face *.node *.halo *.poly *.vtu *.stat *.log-0 *.err-0 */matrixdump*; \
		cd ..; \
	)

