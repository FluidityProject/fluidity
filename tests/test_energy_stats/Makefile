input: clean
	$(call makemesh,channel_2d) ; \
	$(call makemesh,channel_2d_finemesh)

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.face *.halo *.vtu *pvtu *.stat *.log-* *.err-* matrixdump*
	rm -rf lock_exchange_2d_parallel_?
	rm -rf *flredecomp*

makemesh = echo $(1) ; \
           gmsh -2 ./src/$(1).geo ; \
           cp src/$(1).msh .
