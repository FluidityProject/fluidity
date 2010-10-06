XFIG_IMAGES=./numerical_discretisation_images/bisection_method \
./configuring_fluidity_images/periodic_domain \
./numerical_discretisation_images/P1cgshapefunction1d \
./numerical_discretisation_images/P2cgshapefunction1d \
./numerical_discretisation_images/P1cgshapefunction2d \
./numerical_discretisation_images/P2cgshapefunction2d \
./numerical_discretisation_images/P1dgshapefunction1d \
./numerical_discretisation_images/P2dgshapefunction1d \
./numerical_discretisation_images/P1dgshapefunction2d \
./numerical_discretisation_images/P0shapefunction1d \
./numerical_discretisation_images/P0shapefunction2d \
./numerical_discretisation_images/P1cvshapefunction1d \
./numerical_discretisation_images/P2cvshapefunction1d \
./numerical_discretisation_images/P1cvshapefunction2d \
./numerical_discretisation_images/P2cvshapefunction2d \
./numerical_discretisation_images/corner_unstructured_cv \
./numerical_discretisation_images/corner_unstructured \
./numerical_discretisation_images/P1controlvolume3d \
./numerical_discretisation_images/upwind_node_internal \
./numerical_discretisation_images/upwind_node_boundary \
./numerical_discretisation_images/swebynvd \
./numerical_discretisation_images/ultimatenvd \
./numerical_discretisation_images/hypercnvd \
./numerical_discretisation_images/ultracnvd \
./numerical_discretisation_images/ultracmodnvd \
./numerical_discretisation_images/P2dgshapefunction2d

IMAGES=geometry_dimension mesh_name \
./examples_images/lock_exchange/le_basic_0_mesh_nice \
./examples_images/lock_exchange/le_basic_0_T \
./examples_images/lock_exchange/le_basic_10_mesh \
./examples_images/lock_exchange/le_basic_10_T \
./examples_images/lock_exchange/le_basic_30_mesh \
./examples_images/lock_exchange/le_basic_30_T \
./numerical_discretisation_images/bisection_method \
./visualisation_and_diagnostics_images/Paraview-opacity \
./visualisation_and_diagnostics_images/Paraviewconnectbutton \
./visualisation_and_diagnostics_images/Paraviewchooseserver \
./visualisation_and_diagnostics_images/Paraviewchooseserver2 \
./visualisation_and_diagnostics_images/Paraviewconfigureserverwindow2 \
./visualisation_and_diagnostics_images/Paraviewconfigureserverwindow \
./visualisation_and_diagnostics_images/Paraview-filedialogue-timeseries \
./visualisation_and_diagnostics_images/Paraview-rescaling-view \
./visualisation_and_diagnostics_images/Paraview-subwindows \
./visualisation_and_diagnostics_images/diamond_enable_stat  

IMAGE_DIRS=examples_images numerical_discretisation_images \
visualisation_and_diagnostics_images configuring_fluidity_images \
examples_images/sediment\
examples_images/lock_exchange

EPS_IMAGES=./examples_images/rotating_channel/analytic_solution \
./examples_images/rotating_channel/convergence 

CHAPTERS=getting_started model_physics numerical_discretisation parameterisations \
adaptivity mesh_generation configuring_fluidity external_libraries trouble_shooting \
parallel mathematical_notation about visualisation_and_diagnostics examples \
embedded_models paraviewman python dimensionless useful_numbers

fluidity_manual.pdf: fluidity_manual.tex notation.tex\
	$(addsuffix .tex, $(CHAPTERS)) \
	$(addsuffix .pdftex_t, $(XFIG_IMAGES)) \
	$(addsuffix _tex.pdf, $(XFIG_IMAGES))  \
	$(addsuffix .pdf, $(IMAGES))  \
	$(addsuffix .eps, $(EPS_IMAGES))  \
	bibliography.bib

fluidity_manual.html: fluidity_manual.tex notation.tex\
	$(addsuffix .tex, $(CHAPTERS)) \
	$(addsuffix .png, $(XFIG_IMAGES)) \
	$(addsuffix .png, $(IMAGES)) \
	bibliography.bib fluidity_manual.cfg
	htlatex $< "fluidity_manual,2" "" "-dhtml/"
	makeindex -s index.ist $(basename $<)
	if fgrep "Rerun to" $(basename $<.log); then $(MAKE) --assume-new $< $@;fi
ifndef NOBIBTEX
	if fgrep "There were undefined"  $(basename $<.log);\
	then bibtex $(basename $<); \
	$(MAKE) NOBIBTEX=true --assume-new $< $@;fi
endif


#linear_1d_mesh.png: linear_1d_mesh.fig

%.tex: %.pstex_t %.pstex
	../femtools/doc/scripts/wrap_pstex $<

%.dvi: %.tex 
	pdflatex -output-directory $(dir $@) -shell-escape -output-format dvi $<

%.png: %.dvi
	dvipng -T tight  -D 100 -bg Transparent $<  -o $@

PSFILES=$(addsuffix .pdftex_t, $(XFIG_IMAGES)) \
	$(addsuffix _tex.pdf, $(XFIG_IMAGES)) \
	$(addsuffix .png, $(XFIG_IMAGES))\
	$(addsuffix /*.tex, $(IMAGE_DIRS)) \
	$(addsuffix /*.pstex, $(IMAGE_DIRS)) \
	$(addsuffix /*.pstex_t, $(IMAGE_DIRS)) \
	$(addsuffix /*.log, $(IMAGE_DIRS)) \
	$(addsuffix /*.dvi, $(IMAGE_DIRS)) \
	$(addsuffix /*.aux, $(IMAGE_DIRS))

%.ps: %.dvi
	dvips -o $@ $^ 

%.pdf: %.tex
	pdflatex -shell-escape $<
	makeindex -s index.ist $(basename $<)
	if fgrep "Rerun to" $*.log; then $(MAKE) --assume-new $< $@;fi
ifndef NOBIBTEX
	if fgrep "There were undefined"  $*.log;\
	then bibtex $*; \
	$(MAKE) NOBIBTEX=true --assume-new $< $@;fi
endif


%.pdftex: %.fig
	fig2dev -L pdftex $^ $@

%_tex.pdf: %.fig
	fig2dev -L pdftex $^ $@

%.pdftex_t: %.fig
	fig2dev -L pdftex_t -p $*_tex.pdf $^ $@

%.pstex: %.fig
	fig2dev -L pstex $^ $@

%.pstex_t: %.fig
	fig2dev -L pstex_t -p $*.pstex $^ $@

%.pdf: %.fig
	fig2dev -L pdf $^ $@

#%.png: %.fig
#	fig2dev -L png $^ $@

.PHONY: clean

# We don't rm *.ps because that might clobber image files. Instead, put any
# ps files which should go in $(PSFILES)
clean: htclean_fluidity_manual
	rm -r $(TARGETS) $(PSFILES) *.dvi *.aux *.log *.bbl *.blg *.toc *.lof \
*.idx *.ind *.ilg *.out *.pdftex *_tex.pdf *.pdftex_t *.pstex *.pstex_t \
fluidity_manual*.png cmsy*.png $(addsuffix .png, $(XFIG_IMAGES)) 2>/dev/null||true

htclean_%:
	rm -r $*.4ct $*.4tc $*.css $*.idv $*.lg $*.tmp $*.xref *.html 2>/dev/null||true
