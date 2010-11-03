preprocess:
	cp src/circle.geo .
	gmsh -2 circle.geo 
	../../scripts/gmsh2triangle --2d circle.msh; 
	../../bin/fldecomp -n 8  circle;

run:
	@echo **********Calling fluidity with verbose log output enabled:
	../../bin/fluidity -v2 -l rousset-prognostic-new.flml
        
postprocess:
	@echo **********No postprocessing needed

input: preprocess

clean:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *.vtu *.stat *.log-0 *.err-0 */matrixdump*
