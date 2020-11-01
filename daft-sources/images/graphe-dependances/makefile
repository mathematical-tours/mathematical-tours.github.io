# makefile to compile graphe-dependances.png

DOT 		= dot
# FILENAME 	= graphe-codes-correcteurs
# FILENAME 	= graphe-chapitres
FILENAME 	= graphe-dependances

all:	$(FILENAME).ps $(FILENAME).png
ps:	$(FILENAME).ps 
png:	$(FILENAME).png
jpg:	$(FILENAME).jpg
gif:	$(FILENAME).gif

$(FILENAME).ps:			$(FILENAME).dot
	@echo ---
	@echo Compiling ps graph.
	$(DOT) -Tps -o $(FILENAME).eps $(FILENAME).dot

$(FILENAME).png:		$(FILENAME).dot
	@echo ---
	@echo Compiling png graph.
	$(DOT) $(FILENAME).dot -Tpng -o $(FILENAME).png

$(FILENAME).jpg:		$(FILENAME).dot
	@echo ---
	@echo Compiling jpg graph.
	$(DOT) $(FILENAME).dot -Tjpeg -o $(FILENAME).jpg

$(FILENAME).gif:		$(FILENAME).dot
	@echo ---
	@echo Compiling gif graph.
	$(DOT) $(FILENAME).dot -Tgif -o $(FILENAME).gif

clean:

cleanall:
	@echo ---
	@echo Cleaning all.
	rm -f *.eps *.png *.jpg
   
      



