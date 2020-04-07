input:
	N=1; for i in A B C; do \
	  python icosahedral_gmsh.py icosahedron_$${i} 2.22 $${N}; N=$$(($${N}+1)); \
	  sed "s/_TEMPLATE/_$${i}/" spherical_template.flml > spherical_$${i}.flml; \
	done
	spud-set spherical_A.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.25
	spud-set spherical_B.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.125
	spud-set spherical_C.flml "/geometry/mesh::ExtrudedMesh/from_mesh/extrude/layer::WholeDepth/regions::WholeMesh/sizing_function/constant" 0.0625

clean:
	rm -rf spherical_[A-Z].*
	rm -rf fluidity.* *.msh *.vtu *.stat
	rm -rf matrixdump*
	rm -rf *pycache*
