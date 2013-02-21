input: clean
	../../bin/interval --dx 0.05 --reverse -- -0.5 0.5 line
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l mmat-balancepressure-p0p1cv-gpp2 mmat-balancepressure-p0p1cv-gpp2_flredecomp
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l mmat-balancepressure-p0p1cv-hpp2dg mmat-balancepressure-p0p1cv-hpp2dg_flredecomp
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l mmat-balancepressure-p0p1-rhop0-vfracp0 mmat-balancepressure-p0p1-rhop0-vfracp0_flredecomp
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l mmat-balancepressure-p1dgp2-rhop1dg-vfracp1cv mmat-balancepressure-p1dgp2-rhop1dg-vfracp1cv_flredecomp

clean:
	rm -rf *.stat *.vtu *.pvtu *.convergence *.node *.bound *.ele *.halo *.log-? *.err-? *checkpoint*
	rm -rf mmat-balancepressure*_?
	rm -rf mmat-balancepressure*_??
	rm -rf matrixdump*
	rm -rf *flredecomp*
