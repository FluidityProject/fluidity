input: clean
	$(call makemesh,channel_2d) ; \
        ../../bin/fldecomp -n 3 channel_2d
	$(call makemesh,channel_2d_finemesh)

clean: 
	rm -f  ./src/*.msh *.ele *.edge *.node *.face *.halo *.vtu *pvtu *.stat *.log-* *.err-* matrixdump*

makemesh = echo $(1) ; \
           gmsh -2 ./src/$(1).geo ; \
           ../../scripts/gmsh2triangle --2d ./src/$(1).msh
