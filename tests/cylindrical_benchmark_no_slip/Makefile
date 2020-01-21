input:
	N=1; for i in A B C; do \
	  gmsh -2 -o circle_$${i}.msh -setnumber N $${N} mesh/circle.geo; N=$$(($$N*2)); \
	  sed "s/_TEMPLATE/_$${i}/" cylindrical_template.flml > cylindrical_$${i}.flml; \
	done
	spud-set cylindrical_A.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::UpperMantle/regions::WholeMesh/sizing_function/constant" 0.0625
	spud-set cylindrical_A.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::LowerMantle/regions::WholeMesh/sizing_function/constant" 0.0625
	spud-set cylindrical_B.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::UpperMantle/regions::WholeMesh/sizing_function/constant" 0.03125
	spud-set cylindrical_B.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::LowerMantle/regions::WholeMesh/sizing_function/constant" 0.03125
	spud-set cylindrical_C.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::UpperMantle/regions::WholeMesh/sizing_function/constant" 0.015625
	spud-set cylindrical_C.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::LowerMantle/regions::WholeMesh/sizing_function/constant" 0.015625

clean:
	rm -rf cylindrical_[A-Z].*
	rm -rf fluidity.* *.msh *.vtu
	rm -rf matrixdump*
	rm -rf *pycache*
