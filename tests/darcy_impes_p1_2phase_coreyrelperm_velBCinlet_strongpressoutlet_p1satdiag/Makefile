include ../darcy_impes_common/Depend.mk

input: clean
	cp $(FLUIDITYPATH)/tests/darcy_impes_common/mesh_data/line_*.msh .
	cp $(FLUIDITYPATH)/tests/darcy_impes_common/mesh_data/square_*.msh .
	cp $(FLUIDITYPATH)/tests/darcy_impes_common/mesh_data/cube_*.msh .

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
	rm -f error_*.txt
