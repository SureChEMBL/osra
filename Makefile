######################### GCC Setup ###################
X11LIBS=-L/usr/X11R6/lib
#X11LIBS=-L/usr/X11R6/lib64
IMLIBS= $(X11LIBS) -lMagick++ -lWand -lMagick -ltiff -lfreetype -ljpeg -lXext -lSM -lICE -lX11 -lXt -lbz2
#IMLIBS=-L/usr/local/lib -lMagick++ -lWand -lMagick -lgdi32 -ljbig -llcms -ltiff -ljasper  -ljpeg  -lpng -lbz2 -lz 
POTRACELIB=-L../../potrace-1.7/src/ -lpotrace
POTRACEINC=-I../../potrace-1.7/src/
GOCRSRC=../../gocr-0.45/src/
GOCRINC= -I$(GOCRSRC) -I../../gocr-0.45/include/
GOCRLIB= -L$(GOCRSRC) -lPgm2asc -lnetpbm
#OPENBABELLIB=../openbabel-2.1.1/src/formats/cansmilesformat.o ../openbabel-2.1.1/src/formats/obmolecformat.o -L/usr/local/lib -lopenbabel
OPENBABELLIB=-L/usr/local/lib -lopenbabel   
OPENBABELINC=-I/usr/local/include/openbabel-2.0/
TCLAPINC=-I/usr/local/include/tclap/ -I/usr/local/include

OCRAD=../../ocrad-0.16/
OCRADSRC=$(wildcard $(OCRAD)*.cc)
OCRADINC=$(wildcard $(OCRAD)*.h)
OCRADOBJ=$(OCRADSRC:.cc=.o)

CPP = g++ -g -O2 -fPIC -I$(OCRAD) -I/usr/local/include -D_LIB -D_MT -Wall -DHAVE_CONFIG_H $(POTRACEINC) $(GOCRINC) $(OPENBABELINC) $(TCLAPINC)
LD=g++ -g -O2 -fPIC
CP=cp
SED=sed

LIBS=$(POTRACELIB) -lm  $(IMLIBS) $(GOCRLIB) $(OPENBABELLIB) -lz
OBJ = osra.o osra_mol.o  osra_anisotropic.o osra_ocr.o $(OCRADOBJ)


all:	$(OBJ)
	${LD}  -o osra $(OPTS) $(OBJ) $(LIBS)


osra.o:	osra.cpp osra.h pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) -c osra.cpp

osra_mol.o: osra_mol.cpp osra.h
	    $(CPP) -c osra_mol.cpp

osra_anisotropic.o:	osra_anisotropic.cpp osra.h CImg.h greycstoration.h
	$(CPP) -c osra_anisotropic.cpp

osra_ocr.o:	osra_ocr.cpp osra.h $(OCRADSRC) $(OCRADINC) pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) -c osra_ocr.cpp

clean:	
	rm -f *.o osra pgm2asc.h output.h list.h unicode.h gocr.h

pgm2asc.h: $(GOCRSRC)/pgm2asc.h
	$(CP) $(GOCRSRC)/pgm2asc.h ./
output.h: $(GOCRSRC)/output.h
	$(CP) $(GOCRSRC)/output.h ./
gocr.h: $(GOCRSRC)/gocr.h
	$(CP) $(GOCRSRC)/gocr.h ./	
unicode.h: $(GOCRSRC)/unicode.h
	$(SED) '/INFINITY/d' $(GOCRSRC)/unicode.h >unicode.h
list.h: $(GOCRSRC)/list.h
	$(SED) 's/struct\ list/struct\ list\_s/' $(GOCRSRC)/list.h >list.h


$(OCRADOBJ):    $(OCRAD)Makefile $(OCRADSRC) $(OCRADINC) 
	make -C $(OCRAD)   

$(OCRAD)Makefile:  ocrad.Makefile.patch
	cd $(OCRAD);./configure  
	-patch -u -t -N  $(OCRAD)Makefile ocrad.Makefile.patch

$(OCRAD)main.cc: main.cc.patch  
	-patch -u -t -N  $(OCRAD)main.cc main.cc.patch

$(OCRAD)character.h: character.h.patch
	-patch -u -t -N  $(OCRAD)character.h character.h.patch
