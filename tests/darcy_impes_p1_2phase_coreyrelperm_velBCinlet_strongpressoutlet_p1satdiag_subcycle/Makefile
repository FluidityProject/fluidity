SHELL = sh

input: clean
	cp mesh_data/square_A.msh .
	cp ../darcy_impes_common/darcy_impes_tools.py .

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
	rm -f *checkpoint*
	rm -f darcy_impes_tools.py
