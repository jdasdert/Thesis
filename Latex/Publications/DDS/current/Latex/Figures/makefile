SRC = $(wildcard ./*.tex)
PDF = $(SRC:.tex=.pdf)

all: $(SRC)
	make $(^:.tex=.pdf)
	make clean

rDDS_principle.tex: sub/DDS_lowM.pdf sub/DDS_highM.pdf sub/DDS_lowC.pdf sub/DDS_highC.pdf
meshing.pdf: sub/sample_top.pdf
expSim.pdf: sub/samples_small.png sub/longtermSmall_plot.tex sub/comparisonPorosity_plot.tex sub/JPorosity_plot.tex sub/comparisonConcentration_plot.tex sub/JConcentration_plot.tex
influence1.pdf:  sub/releaseRate_plot.tex sub/releaseTime_plot.tex sub/porosity1.pdf sub/Porosity_plot.tex sub/porosity1.pdf sub/porosity2.pdf sub/membraneThickness_plot.tex sub/DDS_thin.pdf sub/DDS_highC.pdf sub/DDS_lowC.pdf sub/DDS_highM.pdf sub/DDS_lowV.pdf sub/concentration_plot.tex sub/volumeReservoir_plot.tex sub/membraneHeight_plot.tex sub/reservoirRadius_plot.tex sub/releaseRate_plot.tex sub/releaseTime_plot.tex
influence2.pdf: sub/CTshift.png sub/centerShift_plot.tex sub/CTcone1.png sub/CTcone2.png sub/tpod_fat.png sub/tpod_slim.png sub/coneShape_plot.tex sub/lrratio_plot.tex sub/sizeTpod_plot.tex
longterm.pdf: sub/longtermRelease_plot.tex sub/samples.png
hydrogel.tex: data/hydrogel.csv data/hydrogel_rates.csv
flowRelease.tex: sub/flowRelease_plot.tex
graphicalAbstract.pdf: abstractExperiments.tex abstractSimulations.pdf abstractParameters.pdf abstractCustomRelease.pdf

%.pdf: %.tex sub/stdSettings.tex makefile
	pdflatex -halt-on-error $<
	pdftoppm -singlefile -png -r 300 $(<:.tex=.pdf) $(<:.tex=)

clean:
	@rm -f *.aux *.log