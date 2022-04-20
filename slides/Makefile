JobName=slides-nnpde-project
TeX=pdflatex
LatexFlags= -jobname=$(JobName) --shell-escape

default:
	make build
	make build
	make clean
	cp $(JobName).pdf ../

build:
	$(TeX) $(LatexFlags) main.tex

clean:
	rm -f $(JobName).log $(JobName).aux $(JobName).out $(JobName).toc *.cpt
	rm -f $(JobName).nav $(JobName).snm $(JobName).vrb
	rm -rf _minted-*
