
DOCS = two-phase-manual.pdf

.PHONY: all

all: $(DOCS)

%.pdf: %.tex 
	pdflatex $< 
	pdflatex $<

