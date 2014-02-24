SHELL = sh

input: clean
	cp ../darcy_impes_common/darcy_impes_tools.py .
	cp ../darcy_impes_common/mesh_data/line_mesh_1d_A* .
	cp ../darcy_impes_common/mesh_data/line_mesh_1d_B* .
	cp ../darcy_impes_common/mesh_data/line_mesh_1d_C* .
	cp ../darcy_impes_common/mesh_data/line_mesh_1d_D* .
	cp ../darcy_impes_common/mesh_data/square_A.msh .
	cp ../darcy_impes_common/mesh_data/square_B.msh .
	cp ../darcy_impes_common/mesh_data/square_C.msh .
	cp ../darcy_impes_common/mesh_data/cube_A.msh .
	cp ../darcy_impes_common/mesh_data/cube_B.msh .

clean:
	rm -f *.ele
	rm -f *.node
	rm -f *.bound
	rm -f *.edge
	rm -f *.face
	rm -f *.vtu 
	rm -f *.pvtu 
	rm -f *.s 
	rm -f *.stat 
	rm -f *.log-0 
	rm -f *.err-0 
	rm -f *.msh
	rm -f *.halo
	rm -f *.err*
	rm -f *.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
	rm -f *.convergence 
	rm -f *.pyc
	rm -f *.txt
	rm -f darcy_impes_tools.py
