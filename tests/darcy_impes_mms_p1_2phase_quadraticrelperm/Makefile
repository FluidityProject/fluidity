include ../darcy_impes_common/Depend.mk

input: code geo clean-generated

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

gen:
	python processing.py gen

xml:
	python processing.py xml

pre:
	python processing.py pre

proc:
	python processing.py proc

post:
	python processing.py post

all:
	python processing.py all

clean-generated:
	python processing.py clean
	rm -f *.vtu 
	rm -f *.pvtu 
	rm -f *.s 
	rm -f *.log-0 
	rm -f *.err-0 
	rm -f *.halo
	rm -f *.err*
	rm -f *.log*
	rm -f matrixdump*
	rm -f first_timestep_adapted_mesh*
	rm -f *.convergence 
	rm -f *.pyc
	rm -f *.stat
	rm -f *.txt

clean-geo:
	rm -f template_*.sh
	rm -f template_*.geo

clean-code:
	rm -f manufactured_solution_test_tools.py
	rm -f test_tools.py
	rm -f batch_tools.py

clean: clean-generated clean-geo clean-code

