DIR = $(wildcard ./*/)
SRC = $(wildcard ./*.tex)
PDF = $(SRC:.tex=.pdf)
PDF_EXTRA = $(filter-out $(PDF),$(wildcard ./*.pdf))
PNG_EXTRA = $(wildcard ./*.png)

all: $(DIR) $(PDF) $(PDF_EXTRA) $(PNG_EXTRA) #$(EPS)

$(DIR): FORCE
	@make -C $@

%.pdf: %.tex
	@pdflatex -halt-on-error $<

.PRECIOUS: $(PDF)

FORCE: ;

clean:
	@rm -f *.dvi *.ps *.aux *.log

allclean:
	@rm -f $(EPS) $(PDF) *-eps-converted-to.pdf *.dvi *.ps *.aux *.log