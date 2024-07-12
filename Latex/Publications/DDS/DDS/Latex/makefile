FIGDIR = ./Figures

pdflatex:
	make figures
	pdflatex -halt-on-error -interaction=batchmode diffModelling.tex
	pdflatex -halt-on-error -interaction=batchmode diffModelling.tex
	#@rm -f *.aux *.log

figures:
	cd $(FIGDIR) && make

clean:
	@rm -f *.aux *.log *.pdf
	cd $(FIGDIR) && make clean