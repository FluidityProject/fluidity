#
input: clean
	cp src/mesh.zip .
	unzip mesh.zip
	triangle2gmsh first_timestep_adapted_mesh

clean: 
	find . -name '*vtu' | xargs -n 10 rm -f
	rm -f *.msh *.halo *.ele *.face *.node *.poly *.s mesh.zip
	rm -rf *checkpoint* 
	rm -f *.stat *.log-* *.err-* matrixdump* *~
	rm -rf lock_exchange_3d_*[0-9]
	rm -rf *flredecomp*
