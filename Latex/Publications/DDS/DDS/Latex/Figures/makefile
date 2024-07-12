SRC = $(wildcard ./*.tex)
PDF = $(SRC:.tex=.pdf)
AX=d025c1000 d030c0050 d045c0025 d045c0050 d045c0100 d045c0200 d050c0050 d055c0050 d060c0050 d070c0050

all: $(SRC)
	make $(^:.tex=.pdf)
	make clean

coneAx.pdf: CTcone1.png CTcone2.png
cshiftAx.pdf: CTshift.png
expComparison.pdf: $(AX:=Ax.pdf)
compConc.pdf: Jconc.pdf comparisonConc.pdf
compDens.pdf: Jdens.pdf comparisonDens.pdf
compDensMB_AB.pdf: JdensMB_AB.pdf comparisonDensMB_AB.pdf
simComparison.pdf: ZnOdensAx.pdf thicknessAx.pdf cshiftAx.pdf coneAx.pdf lrratioAx.pdf volumeTpodAx.pdf volumeResAx.pdf concentrationAx.pdf

%Ax.pdf: %Ax.tex %.pdf makefile
	pdflatex -halt-on-error -interaction=batchmode $<
	pdftoppm -singlefile -png -r 1200 $(<:.tex=.pdf) $(<:.tex=)

%.pdf: %.tex makefile
	pdflatex -halt-on-error -interaction=batchmode $<
	pdftoppm -singlefile -png -r 1200 $(<:.tex=.pdf) $(<:.tex=)

clean:
	@rm -f *.aux *.log

veryclean:
	make clean
	@rm -f $(SRC:.tex=.pdf)