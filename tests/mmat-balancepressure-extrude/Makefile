input: clean
	../../bin/interval --dx 0.05 --reverse -- -0.5 0.5 line
	../../bin/fldecomp -n 2 line

clean:
	rm -rf *.halo *.log *.d.* *.stat *.vtu *.pvtu *.convergence *.node *.bound *.ele *checkpoint*
