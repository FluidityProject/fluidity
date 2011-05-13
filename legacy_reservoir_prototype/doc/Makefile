
default:	Model_Documentation.pdf 

Model_Documentation.pdf:	Model_Documentation.tex
Model_Documentation.idx:	Model_Documentation.idx

%.pdf:	%.tex
	rm -f $*.pdf
	pdflatex $<
	makeindex $*.idx
	pdflatex $<
	if fgrep "Rerun to" $*.log; then $(MAKE) --assume-new $< $@;fi
ifndef NOBIBTEX
	if fgrep "There were undefined"  $*.log;\
	then bibtex $*; \
	$(MAKE) NOBIBTEX=true --assume-new $< $@;fi
endif

%.idx:	%.idx
	makeindex $<

.PHONY: clean clean-pdf

# We don't rm *.ps because that might clobber image files. Instead, put any
# ps files which should go in $(PSFILES)
clean: htclean_QCCSRC
	rm $(TARGETS)*.aux *.log *.bbl *.blg *.toc *.lof *ind *idx *ilg \
*.out *.pdftex *_tex.pdf *.nav *.snm *.lot *.tex~ *.bib~ \
*.pdftex_t *.pstex *.pstex_t Model_Documentation.pdf \
 cmsy*.png $(addsuffix .png, $(XFIG_IMAGES)) 2>/dev/null||true

htclean_%:
	rm $*.4ct $*.4tc $*.css $*.idv $*.lg $*.tmp $*.xref $*.html $*.nav $*.snm \
$*.lot $*.tex~ $*.bib~ 2>/dev/null||true

clean-pdf: 
		rm Model_Documentation.pdf 
