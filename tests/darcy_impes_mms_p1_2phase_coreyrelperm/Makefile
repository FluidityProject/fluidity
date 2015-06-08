include ../darcy_impes_common/Depend.mk

all: geo
	python processing.py all

input: geo clean-generated

geo:
	cp $(FLUIDITYPATH)/tests/darcy_impes_common/mesh_data/template_* .

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

clean: clean-generated clean-geo

