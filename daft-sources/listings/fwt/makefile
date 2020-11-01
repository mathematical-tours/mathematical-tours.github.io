# makfile pour compiler les listings de FWT.
ML2LATEX = perl ../ml2latex.pl 
FILE1 = fwt

all: $(FILE1).tex

$(FILE1).tex:	$(FILE1).m
	$(ML2LATEX)	$(FILE1).m