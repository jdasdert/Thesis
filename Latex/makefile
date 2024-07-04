MAIN = Dittmann_Thesis.tex
LFLAGS = -halt-on-error -interaction=nonstopmode

main: 
	@make pdflatex

pdflatex: $(MAIN)
	@make -C Figs
	@pdflatex $(LFLAGS) $(MAIN) 
	@biber $(MAIN:.tex=)
	@pdflatex $(LFLAGS) $(MAIN)
	@pdflatex $(LFLAGS) $(MAIN)
	@make clean

latex: $(MAIN)
	@make -C Figs
	@latex $(MAIN)

clean:
	@make -C Figs clean
	@rm -f *.aux *.bbl *.blg *.dvi *.log *.ps

allclean:
	@make -C Figs allclean
	@rm -f *.pdf *.log *.aux