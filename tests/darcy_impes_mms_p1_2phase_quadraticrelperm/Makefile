include ../darcy_impes_common/Depend.mk

default: clean code geo

code:
	$(call import_python_code,batch_tools.py)
	$(call import_python_code,test_tools.py)
	$(call import_python_code,manufactured_solution_test_tools.py)

geo:
	cp ../darcy_impes_common/mesh_data/template_line.sh .
	cp ../darcy_impes_common/mesh_data/template_rectangle_reg.geo .
	cp ../darcy_impes_common/mesh_data/template_rectangle_irreg.geo .
	cp ../darcy_impes_common/mesh_data/template_cuboid_reg.geo .
	cp ../darcy_impes_common/mesh_data/template_cuboid_irreg.geo .

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
	rm -f *.geo
	rm -f *.sh
	rm -f *.halo
	rm -f *.err*
	rm -f *.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
	rm -f *.convergence 
	rm -f *.pyc
	rm -f error_*.txt
	rm -f manufactured_solution_test_tools.py
	rm -f test_tools.py
	rm -f batch_tools.py
