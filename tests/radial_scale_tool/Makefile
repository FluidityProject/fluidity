input: clean
	@echo **********Converting the gmsh mesh to triangle format:
	../../bin/gmsh2triangle spherical_patch.msh; 
	@echo **********Applying sed script:
	sed '1s/\<3\>/2/' < spherical_patch.node > mod_spherical_patch.node
	mv mod_spherical_patch.node spherical_patch.node

#The above sed script lets fluidity know that this is
#a 2d manifold in 3d space.

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -f *.halo *.vtu *.stat matrixdump* *log* *err* *.node *.ele *.face *.edge *.geo