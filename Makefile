ARCH=unix
#unix,win32

OPENBABEL_OR_RDKIT=rdkit
#openbabel,rdkit

POTRACE=../potrace-1.8/
GOCR=../gocr-0.45/
OCRAD=../ocrad-0.18-pre4/

ifeq ($(OPENBABEL_OR_RDKIT),rdkit)
#RDKIT=../RDKit_Q32008_1/
RDKIT=../rdkit-svn/
BOOST=../boost_1_37_0/
endif

ifeq ($(ARCH),win32)
OPENBABEL=../openbabel-2.2.0/
endif


CPP = g++
LD=g++ -g -O2 -fPIC
CP=cp
SED=sed
RM=rm
PATCH=patch

SRCDIR=./

################ Hopefully you won't have to change anything below this line #########

TCLAPINC=-I/usr/local/include/tclap/ -I/usr/local/include

RDKITINC=-I$(RDKIT)/Code/ -I$(RDKIT)/External/vflib-2.0/include/ -I$(BOOST)
RDKITLIB=-L$(RDKIT)/bin/ -lRDGeneral -lSmilesParse -lGraphMol -lFileParsers -lDepictor -lRDGeometry -lSubstruct -lDistGeomHelpers -lDistGeom -lForceFieldHelpers -lForceField -lEigenSolvers -lAlignment -lOptimizer -L$(RDKIT)/External/vflib-2.0/lib -lvf

ifneq ($(ARCH),win32)
NETPBM=-lnetpbm
else
MAGIKLIB_WIN32=-lgdi32 -llcms -ljbig  -ltiff -ljasper  -ljpeg -lpng
endif

MAGIKINC := $(shell Magick++-config --cppflags) $(shell Magick++-config --cxxflags)
MAGIKLIB := $(shell Magick++-config --libs) $(shell Magick++-config --ldflags) $(MAGIKLIB_WIN32)


POTRACELIB=-L$(POTRACE)/src/ -lpotrace
POTRACEINC=-I$(POTRACE)/src/
GOCRSRC=$(GOCR)/src/
GOCRINC= -I$(GOCRSRC) -I$(GOCR)/include/
GOCRLIB= -L$(GOCRSRC) -lPgm2asc $(NETPBM)

ifeq ($(ARCH),win32)
#OPENBABELLIB=$(OPENBABEL)/src/formats/cansmilesformat.o $(OPENBABEL)/src/formats/obmolecformat.o -L/usr/local/lib -lopenbabel
OPENBABELLIB=-L/usr/local/lib -lopenbabel
OPENBABELINC=-I$(OPENBABEL)/include
else
OPENBABELLIB=-L/usr/local/lib -lopenbabel   
OPENBABELINC=-I/usr/include/openbabel-2.0/ -I/usr/local/include/openbabel-2.0/
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


OCRADSRC=$(wildcard $(OCRAD)*.cc)
OCRADINC=$(wildcard $(OCRAD)*.h)
OCRADOBJ=$(OCRADSRC:.cc=.o)

CPPFLAGS= -g -O2 -fPIC -I$(OCRAD) -I/usr/local/include -D_LIB -D_MT -Wall $(POTRACEINC) $(GOCRINC) $(MOL_BACKEND_INC) $(TCLAPINC) $(MAGIKINC)

LIBS=$(POTRACELIB) -lm  $(MAGIKLIB) $(GOCRLIB) $(MOL_BACKEND_LIB) -lz
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
