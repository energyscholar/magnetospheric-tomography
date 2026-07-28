TEX = Stephenson_AltitudeResolvedCoherence_2026
PDF = $(TEX).pdf

all: $(PDF)

$(PDF): $(TEX).tex $(TEX).bib figures/*.pdf
	pdflatex $(TEX)
	bibtex $(TEX)
	pdflatex $(TEX)
	pdflatex $(TEX)

clean:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out $(TEX).pdf

.PHONY: all clean
