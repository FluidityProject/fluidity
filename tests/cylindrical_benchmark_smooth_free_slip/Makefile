input:
	N=1; for i in A B C D E; do \
	  gmsh -2 -o circle_$${i}.msh -setnumber N $${N} src/circle.geo; N=$$(($$N*2)); \
	  sed "s/_TEMPLATE/_$${i}/" cylindrical_template.flml > cylindrical_$${i}.flml; \
	done
	spud-set cylindrical_A.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.125
	spud-set cylindrical_B.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.0625
	spud-set cylindrical_C.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.03125
	spud-set cylindrical_D.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.015625
	spud-set cylindrical_E.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.0078125



clean:
	rm -rf cylindrical_[A-Z].*
	rm -rf fluidity.* *.msh *.vtu
	rm -rf matrixdump*
