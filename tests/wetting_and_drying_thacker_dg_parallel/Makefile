all: input
	mpiexec flredecomp -i 1 -o 4 -v -l thacker_B_d0=1.0 thacker_B_d0=1.0_flredecomp
	mpiexec -n 4 fluidity -v5 thacker_B_d0=1.0_flredecomp.flml &> output.txt

input: clean
	#cd mesh; ./makemesh.sh

clean:
	rm -f *.halo *.poly *vtu *.s *.d.1 *.stat *.msh *.detectors matrixdump* *.pdf *.dat *.pyc *.err-* *.log-*
	rm -rf thacker_B_d0=1.0_*
	#cd mesh; rm -f  *.ele *.edge *.node *.poly *.vtu *.s *.d.1 *.stat *.detectors matrixdump* 
	rm -rf *flredecomp*
