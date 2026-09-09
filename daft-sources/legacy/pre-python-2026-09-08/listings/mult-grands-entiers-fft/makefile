# makfile pour compiler les listings de multiplication d'entier par Matlab.
ML2LATEX = perl ../ml2latex.pl 
FILE1 = mult_entiers
FILE2 = number2vector
FILE3 = test_mult_entiers
FILE4 = vector2number

all: $(FILE1).tex $(FILE2).tex $(FILE3).tex $(FILE4).tex

$(FILE1).tex:	$(FILE1).m
	$(ML2LATEX)	$(FILE1).m
$(FILE2).tex:	$(FILE2).m
	$(ML2LATEX)	$(FILE2).m
$(FILE3).tex:	$(FILE3).m
	$(ML2LATEX)	$(FILE3).m
$(FILE4).tex:	$(FILE4).m
	$(ML2LATEX)	$(FILE4).m