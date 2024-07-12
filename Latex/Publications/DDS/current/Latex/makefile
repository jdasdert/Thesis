FIGDIR = ./Figures

pdflatex:
	make figures
	pdflatex -halt-on-error DDS.tex
	pdflatex -halt-on-error DDS.tex
	pdflatex -halt-on-error SI.tex
	pdflatex -halt-on-error SI.tex
	@rm -f *.aux *.log

figures:
	cd $(FIGDIR) && make

clean:
	@rm -f *.aux *.log *.pdf
	cd $(FIGDIR) && make clean