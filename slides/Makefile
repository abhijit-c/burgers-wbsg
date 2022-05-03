JobName=slides-nnpde-project
TeX=pdflatex
LatexFlags= -jobname=$(JobName) --output-directory=tmp --shell-escape

default:
	make build
	biber tmp/$(JobName)
	make build
	#make clean
	cp tmp/$(JobName).pdf .

build:
	$(TeX) $(LatexFlags) main.tex

clean:
	rm -f tmp/*
	rm -rf _minted-*
