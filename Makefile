ARCH=i386
#i386,x86_64,win32,osx

OPENBABEL_OR_RDKIT=openbabel
#openbabel,rdkit

POTRACE=../../potrace-1.7/
GOCR=../../gocr-0.45/
OCRAD=../../ocrad-0.17/
RDKIT=../../RDKit_64/

ifeq ($(ARCH),win32)
OPENBABEL=../openbabel-2.1.1/
endif

ifeq ($(ARCH),osx)
MAGIK=/usr/ImageMagick-6.3.5/
endif


CPP = g++
LD=g++ -g -O2 -fPIC
CP=cp
SED=sed
RM=rm
PATCH=patch

SRCDIR=./

RDKITINC=-I$(RDKIT)/Code/ -I$(RDKIT)/External/vflib-2.0/include/ -I../../boost_1_34_1
RDKITLIB=-L$(RDKIT)/bin/ -lRDGeneral -lSmilesParse -lGraphMol -lFileParsers -lDepictor -lRDGeometry -lSubstruct -L$(RDKIT)/External/vflib-2.0/lib -lvf

################ Hopefully you won't have to change anything below this line #########

ifeq ($(ARCH),x86_64)
X11LIBS=-L/usr/X11R6/lib64
else
X11LIBS=-L/usr/X11R6/lib
endif

ifeq ($(ARCH),osx)
MAGIKINC=-I$(MAGIK)/include/
MAGIKLIB=-L$(MAGIK)/lib/  -lMagick++ -lWand -lMagick -ltiff -ljpeg
IMLIBS= $(X11LIBS) $(MAGIKLIB) -lfreetype -lXext -lSM -lICE -lX11 -lXt -lbz2 -lz
endif

ifeq ($(ARCH),win32)
IMLIBS=-L/usr/local/lib -lMagick++ -lWand -lMagick -lgdi32 -ljbig -llcms -ltiff -ljasper  -ljpeg  -lpng -lbz2 -lz 
else
IMLIBS= $(X11LIBS) -lMagick++ -lWand -lMagick -ltiff -lfreetype -ljpeg -lXext -lSM -lICE -lX11 -lXt -lbz2
NETPBM=-lnetpbm
endif

POTRACELIB=-L$(POTRACE)/src/ -lpotrace
POTRACEINC=-I$(POTRACE)/src/
GOCRSRC=$(GOCR)/src/
GOCRINC= -I$(GOCRSRC) -I$(GOCR)/include/
GOCRLIB= -L$(GOCRSRC) -lPgm2asc $(NETPBM)

ifeq ($(ARCH),win32)
OPENBABELLIB=$(OPENBABEL)/src/formats/cansmilesformat.o $(OPENBABEL)/src/formats/obmolecformat.o -L/usr/local/lib -lopenbabel
OPENBABELINC=-I$(OPENBABEL)/include
else
OPENBABELLIB=-L/usr/local/lib -lopenbabel   
OPENBABELINC=-I/usr/local/include/openbabel-2.0/
endif

ifeq ($(OPENBABEL_OR_RDKIT),rdkit)
MOL_BACKEND_INC=$(RDKITINC)
MOL_BACKEND_LIB=$(RDKITLIB)
MOL_BACKEND_CPP=osra_rdkit.cpp
MOL_BACKEND_OBJ=osra_rdkit.o
LD_LIBRARY_PATH=$(RDKIT)/bin
else
MOL_BACKEND_INC=$(OPENBABELINC)
MOL_BACKEND_LIB=$(OPENBABELLIB)
MOL_BACKEND_CPP=osra_openbabel.cpp
MOL_BACKEND_OBJ=osra_openbabel.o
endif

TCLAPINC=-I/usr/local/include/tclap/ -I/usr/local/include

OCRADSRC=$(wildcard $(OCRAD)*.cc)
OCRADINC=$(wildcard $(OCRAD)*.h)
OCRADOBJ=$(OCRADSRC:.cc=.o)

CPPFLAGS= -g -O2 -fPIC -I$(OCRAD) -I/usr/local/include -D_LIB -D_MT -Wall $(POTRACEINC) $(GOCRINC) $(MOL_BACKEND_INC) $(TCLAPINC) $(MAGIKINC)

LIBS=$(POTRACELIB) -lm  $(IMLIBS) $(GOCRLIB) $(MOL_BACKEND_LIB) -lz
OBJ = osra.o osra_anisotropic.o osra_ocr.o $(MOL_BACKEND_OBJ) $(OCRADOBJ)


all:	$(OBJ)
	${LD}  -o osra $(OBJ) $(LIBS)


osra.o:	osra.cpp osra.h pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) $(CPPFLAGS) -c osra.cpp

$(MOL_BACKEND_OBJ): $(MOL_BACKEND_CPP) osra.h
	$(CPP) $(CPPFLAGS) -c $(MOL_BACKEND_CPP)


osra_anisotropic.o:	osra_anisotropic.cpp osra.h CImg.h greycstoration.h
	$(CPP) $(CPPFLAGS) -c osra_anisotropic.cpp

osra_ocr.o:	osra_ocr.cpp osra.h $(OCRADSRC) $(OCRADINC) pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) $(CPPFLAGS) -c osra_ocr.cpp

clean:	
	-$(RM) -f *.o osra pgm2asc.h output.h list.h unicode.h gocr.h

pgm2asc.h: $(GOCRSRC)/pgm2asc.h
	$(CP) $(GOCRSRC)/pgm2asc.h $(SRCDIR)
output.h: $(GOCRSRC)/output.h
	$(CP) $(GOCRSRC)/output.h $(SRCDIR)
gocr.h: $(GOCRSRC)/gocr.h
	$(CP) $(GOCRSRC)/gocr.h $(SRCDIR)	
unicode.h: $(GOCRSRC)/unicode.h
	$(SED) '/INFINITY/d' $(GOCRSRC)/unicode.h >$(SRCDIR)/unicode.h
list.h: $(GOCRSRC)/list.h
	$(SED) 's/struct\ list/struct\ list\_s/' $(GOCRSRC)/list.h >$(SRCDIR)/list.h


$(OCRADOBJ):    $(OCRAD)/Makefile $(OCRADSRC) $(OCRADINC) 
	$(MAKE) -C $(OCRAD)   

$(OCRAD)/Makefile:  ocrad.Makefile.patch
	cd $(OCRAD);./configure  
	-$(PATCH) -u -t -N  $(OCRAD)/Makefile $(SRCDIR)/ocrad.Makefile.patch

$(OCRAD)/main.cc: main.cc.patch  
	-$(PATCH) -u -t -N  $(OCRAD)/main.cc $(SRCDIR)/main.cc.patch

$(OCRAD)/character.h: character.h.patch
	-$(PATCH) -u -t -N  $(OCRAD)/character.h $(SRCDIR)/character.h.patch

$(OCRADSRC): $(OCRAD)/main.cc

$(OCRADINC): $(OCRAD)/character.h
