ARCH=unix
#unix,win32, osx-static

POTRACE=../potrace-1.8/
GOCR=../gocr-0.45/
OCRAD=../ocrad-0.19-pre1/
OPENBABEL=/usr/
MAGICKCONFIG=GraphicsMagick++-config

TCLAPINC=-I/usr/local/include/tclap/

ifeq ($(ARCH), osx-static)
STATIC_LIB_DIR=/Users/igor/build/lib
endif

#TESSERACTINC=-I/usr/local/include
#TESSERACTLIB=-L/usr/local/lib -ltesseract_full

CPP = g++  -g -O3
LD=g++ -g -O3 -fPIC
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

MAGIKINC := $(shell $(MAGICKCONFIG) --cppflags) $(shell $(MAGICKCONFIG) --cxxflags)
MAGIKLIB := $(shell $(MAGICKCONFIG) --libs) $(shell $(MAGICKCONFIG) --ldflags) $(MAGIKLIB_WIN32)
ifeq ($(ARCH), osx-static)
NETPBM=
LDFLAGS_STATIC=-static-libgcc -static-libstc++
MAGIKINC := $(MAGIKINC) -Dcimg_display_type=0
MAGIKLIB := -L$(STATIC_LIB_DIR)/ -lGraphicsMagick++ -lGraphicsMagick -llcms -ltiff -lfreetype -ljasper -ljpeg -lpng  $(STATIC_LIB_DIR)/libbz2.a $(STATIC_LIB_DIR)/libz.a -lm
endif

POTRACELIB=-L$(POTRACE)/src/ -lpotrace
POTRACEINC=-I$(POTRACE)/src/
GOCRSRC=$(GOCR)/src/
GOCRINC= -I$(GOCRSRC) -I$(GOCR)/include/
GOCRLIB= -L$(GOCRSRC) -lPgm2asc $(NETPBM)
OCRADINC=-I$(OCRAD)
OCRADLIB=-L$(OCRAD) -locrad

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



CPPFLAGS= -g -O2 -fPIC $(OCRADINC) $(MINGWINC) -D_LIB -D_MT -Wall $(POTRACEINC) $(GOCRINC) $(MOL_BACKEND_INC) $(TCLAPINC) $(MAGIKINC) $(TESSERACTINC)

LIBS=$(POTRACELIB) $(OCRADLIB) -lm  $(MAGIKLIB) $(GOCRLIB) $(MOL_BACKEND_LIB)  $(TESSERACTLIB)
OBJ = osra.o osra_anisotropic.o osra_ocr.o $(MOL_BACKEND_OBJ) $(MCDLUTIL) unpaper.o


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

osra_ocr.o:	osra_ocr.cpp osra.h pgm2asc.h output.h list.h unicode.h gocr.h
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


unpaper.o:	unpaper.c
	$(CPP) $(CPPFLAGS) -c unpaper.c