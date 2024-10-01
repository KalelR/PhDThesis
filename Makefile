MAIN=main
OUTPUT=$(MAIN).pdf

main: $(MAIN).tex
	@pdflatex $(MAIN)
	@bibtex   $(MAIN)
	@pdflatex $(MAIN)
	@pdflatex $(MAIN)
	@mv $(MAIN).pdf $(OUTPUT)

bib:
	@bibtex $(MAIN)

# dessa forma somente os arquivos removidos são listados na saída
clean:
	rm paper-malleability/*.aux
	rm *.aux
	rm 'main.synctex(busy)'
open:
	xdg-open $(OUTPUT)&

purge: clean
	@rm -f $(OUTPUT)

