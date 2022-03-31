all: input
	fluidity  -v5 balzano1_2plus1.flml &> output.txt

input: clean
	cd mesh; ./makemesh.sh

clean:
	rm -f  *.poly *vtu *.s *.d.1 *.stat *.msh *.halo *.detectors matrixdump* *.pdf *.png *.dat
	cd mesh; rm -f *.poly *.vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.halo
	rm -rf *flredecomp*
config:
	diamond -s ../../tools/fluidity_options.rng  balzano1_2plus1.flml 

