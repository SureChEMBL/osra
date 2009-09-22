ARCH=unix
#unix,win32, osx-static

POTRACE=../potrace-1.8/
GOCR=../gocr-0.45/
OCRAD=../ocrad-0.18/
OPENBABEL=/usr/local/

ifeq ($(ARCH), osx-static)
MAGICK_STATIC_INC=-I/Users/igor/build/include/ImageMagick
STATIC_LIB_DIR=/Users/igor/build/lib
endif

TCLAPINC=-I/usr/local/include/tclap/

TESSERACTINC=-I/usr/local/include
TESSERACTLIB=-L/usr/local/lib -ltesseract_full

CPP = g++
LD=g++ -g -O2 -fPIC
CP=cp
SED=sed
RM=rm
PATCH=patch

SRCDIR=./


################ Hopefully you won't have to change anything below this line #########


ifneq ($(ARCH),win32)
NETPBM=-lnetpbm
else
MAGIKLIB_WIN32=-lgdi32 -llcms -ljbig  -ltiff -ljasper  -ljpeg -lpng
MINGWINC=-I/usr/local/include
endif

ifneq ($(ARCH), osx-static)
MAGIKINC := $(shell Magick++-config --cppflags) $(shell Magick++-config --cxxflags)
MAGIKLIB := $(shell Magick++-config --libs) $(shell Magick++-config --ldflags) $(MAGIKLIB_WIN32)
else
NETPBM=
LDFLAGS_STATIC=-static-libgcc -static-libstc++
MAGIKINC := $(MAGICK_STATIC_INC) -Dcimg_display_type=0
MAGIKLIB := -L$(STATIC_LIB_DIR) -lMagick++ -lMagickWand -lMagickCore -ltiff -lfreetype -ljasper -ljpeg -lpng $(STATIC_LIB_DIR)/libbz2.a $(STATIC_LIB_DIR)/libz.a -lm -lfreetype
endif

POTRACELIB=-L$(POTRACE)/src/ -lpotrace
POTRACEINC=-I$(POTRACE)/src/
GOCRSRC=$(GOCR)/src/
GOCRINC= -I$(GOCRSRC) -I$(GOCR)/include/
GOCRLIB= -L$(GOCRSRC) -lPgm2asc $(NETPBM)

ifneq ($(ARCH),win32)
OPENBABELLIB=-L$(OPENBABEL)/lib -lopenbabel
OPENBABELINC=-I$(OPENBABEL)/include/openbabel-2.0/
else
OPENBABELLIB=-L/usr/local/lib -lopenbabel
OPENBABELINC=-I$(OPENBABEL)/include/
endif

MOL_BACKEND_INC=$(OPENBABELINC)
MOL_BACKEND_LIB=$(OPENBABELLIB)
MOL_BACKEND_CPP=osra_openbabel.cpp
MOL_BACKEND_OBJ=osra_openbabel.o
MCDLUTIL=mcdlutil.o


OCRADSRC=$(wildcard $(OCRAD)*.cc)
OCRADINC=$(wildcard $(OCRAD)*.h)
OCRADOBJ=$(OCRADSRC:.cc=.o)


CPPFLAGS= -g -O2 -fPIC -I$(OCRAD) $(MINGWINC) -D_LIB -D_MT -Wall $(POTRACEINC) $(GOCRINC) $(MOL_BACKEND_INC) $(TCLAPINC) $(MAGIKINC) $(TESSERACTINC)

LIBS=$(POTRACELIB) -lm  $(MAGIKLIB) $(GOCRLIB) $(MOL_BACKEND_LIB)  $(TESSERACTLIB)
OBJ = osra.o osra_anisotropic.o osra_ocr.o $(MOL_BACKEND_OBJ) $(OCRADOBJ) $(MCDLUTIL)


all:	$(OBJ)
	${LD} $(LDFLAGS_STATIC)  -o osra $(OBJ) $(LIBS)


osra.o:	osra.cpp osra.h pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) $(CPPFLAGS) -c osra.cpp

$(MOL_BACKEND_OBJ): $(MOL_BACKEND_CPP) osra.h
	$(CPP) $(CPPFLAGS) -c $(MOL_BACKEND_CPP)

$(MCDLUTIL):	mcdlutil.cpp mcdlutil.h
	$(CPP) $(CPPFLAGS) -c mcdlutil.cpp

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
