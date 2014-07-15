include ../darcy_impes_common/Depend.mk

input: clean
	$(call import_python_code,batch_tools.py)
	$(call import_python_code,test_tools.py)
	$(call import_python_code,buckley_leverett_test_tools.py)
	cp mesh_data/square_A.msh .

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
	rm -f error_*.txt
	rm -f buckley_leverett_test_tools.py
	rm -f test_tools.py
	rm -f batch_tools.py
