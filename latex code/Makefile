NAME = tesi
TEX  = $(NAME).tex
DVI  = $(TEX:.tex=.dvi)
PS   = $(DVI:.dvi=.ps)
PDF  = $(PS:.ps=.pdf)
#BIB  = biblio
#BIBTEX = $(BIB).tex
#BIBBIB = $(BIB).bib
#BBL = $(NAME:.tex=.bbl)
#CONVERSION = convertbiblio.py

all: $(PDF)

$(PDF): $(TEX)
	pdflatex $(TEX)

$(PS): $(DVI)
	dvips $(DVI)

$(DVI): $(TEX)
	latex $(TEX)

#$(BIB): $(DVI) biblio.tex
#	$(CONVERSION) $(BIBTEX) > $(BIBBIB)
#	bibtex $(NAME)

clean:
	rm *~ pdfa.xmpi *.xmpdata *.toc *.blg *.log *.aux *eps-converted-to.pdf $(NAME).out $(DVI) $(PS) $(BIB)

