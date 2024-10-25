MAIN=main
OUTPUT=$(MAIN).pdf

all: clean main open

main: $(MAIN).tex
	@pdflatex $(MAIN)
	@bibtex   $(MAIN)
	@pdflatex $(MAIN)
	@pdflatex $(MAIN)

bib:
	@bibtex $(MAIN)

# dessa forma somente os arquivos removidos são listados na saída
clean:
	@rm -f paper-malleability/*.aux
	@rm -f paper-metastability/*.aux
	@rm -f *.aux
	@rm -f 'main.synctex(busy)'
	@rm -f main.bbl
open:
	@xdg-open $(OUTPUT)&


