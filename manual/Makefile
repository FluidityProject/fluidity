XFIG_IMAGES=./numerical_discretisation_images/bisection_method \
./examples_images/sediment/bathymetry \
./configuring_fluidity_images/periodic_domain


IMAGES=geometry_dimension mesh_name \
./examples_images/lock_exchange/le_basic_0_mesh_nice \
./examples_images/lock_exchange/le_basic_0_T \
./examples_images/lock_exchange/le_basic_10_mesh \
./examples_images/lock_exchange/le_basic_10_T \
./examples_images/lock_exchange/le_basic_30_mesh \
./examples_images/lock_exchange/le_basic_30_T \
./numerical_discretisation_images/bisection_method \
./examples_images/sediment/bathymetry \
./examples_images/sediment/sed_conc_20 \
./examples_images/sediment/sed_conc_35 \
./examples_images/sediment/sed_conc_50 \
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

EPS_IMAGES=./examples_images/rotating_channel/analytic_solution \
./examples_images/rotating_channel/convergence 

CHAPTERS=getting_started model_physics numerical_discretisation parameterisations \
adaptivity mesh_generation configuring_fluidity external_libraries trouble_shooting \
parallel mathematical_notation about visualisation_and_diagnostics examples \
embedded_models paraviewman



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
	bibliography.bib
	htlatex $< "fluidity_manual"
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
	$(addsuffix .png, $(XFIG_IMAGES)) 


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
