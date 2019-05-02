# Nauty's source code folder
NAUTY=../nauty/

all: 
	gcc -lm -O4 -I$(NAUTY) -march=native gencrn.c $(NAUTY)gtnauty.c $(NAUTY)gtools.c $(NAUTY)naugraph.c $(NAUTY)naugroup.c $(NAUTY)naurng.c $(NAUTY)nausparse.c $(NAUTY)nautil.c $(NAUTY)nautinv.c $(NAUTY)naututil.c $(NAUTY)nauty.c $(NAUTY)schreier.c $(NAUTY)naututil.h -o genCRN
	gcc genInputGraphs.c -o genInputGraphs
