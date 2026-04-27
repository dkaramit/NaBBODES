TEX=pdflatex
BIB=bibtex

FLG=-shell-escape -synctex=1

#make the pdf
all:  Basics.pdf

#run pdflatex and bibtex (runs only if you update something in pdfs or in the tex file)
Basics.pdf:  Basics.tex
	$(TEX) $(FLG) $< 
	$(TEX) $(FLG) $< 
# 	$(BIB) Basics || true
# 	$(TEX) $(FLG) $< 
# 	$(TEX) $(FLG) $< 

	#remove all output except the pdf file
	rm -f *.aux 
	rm -f *.log 
	rm -f *.out 
	rm -f *.bbl 
	rm -f *.toc 
	rm -f *.blg 
	rm -f *.lol
	rm -rf _minted-Basics
clean:
	rm -f Basics.pdf 
	rm -f Basics.pdf
	rm -f *.aux 
	rm -f *.log 
	rm -f *.out 
	rm -f *.bbl 
	rm -f *.blg 
	rm -f *.toc 
	rm -f *.lol
	rm -rf _minted-Basics
	