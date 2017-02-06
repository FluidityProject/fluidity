input:
	N=1; for i in A B C D; do \
	  gmsh -2 -o transfinite_cylinder_$${i}.msh -setnumber N $${N} src/transfinite_cylinder.geo; N=$$(($$N*2)); \
	  sed "s/_TEMPLATE/_$${i}/" cylindrical_template.flml > cylindrical_$${i}.flml; \
	done



clean:
	rm -rf cylindrical_[A-Z].*
	rm -rf fluidity.* *.msh *.vtu
	rm -rf matrixdump*
