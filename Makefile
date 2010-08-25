ARCH=unix
#unix,win32, osx-static
#TESSERACT_ENABLE=yes
CUNEIFORM_ENABLE=yes


POTRACE=../potrace-1.8/
GOCR=../gocr-0.45-patched/
OCRAD=../ocrad-0.20/
OPENBABEL=/usr/
MAGICKCONFIG=GraphicsMagick++-config

TCLAPINC=-I/usr/local/include/tclap/

ifeq ($(ARCH), osx-static)
STATIC_LIB_DIR=/Users/igor/build/lib
endif

ifeq ($(TESSERACT_ENABLE),yes)
TESSERACTINC=-I/usr/local/include
TESSERACTLIB=-L/usr/local/lib -ltesseract_full -ltiff
endif

ifeq ($(CUNEIFORM_ENABLE),yes)
CUNEIFORMINC=-I../cuneiform-linux-1.0.0/install/include -I../cuneiform-linux-1.0.0/cuneiform_src/Kern/hhh/ -I../cuneiform-linux-1.0.0/cuneiform_src/Kern/h/  -I../cuneiform-linux-1.0.0/cuneiform_src/Kern/hrk/  -I../cuneiform-linux-1.0.0/cuneiform_src/Kern/puma/h/    -I../cuneiform-linux-1.0.0/cuneiform_src/Kern/include/ -I../cuneiform-linux-1.0.0/cuneiform_src/Kern/hdebug/
CUNEIFORMLIB=-L../cuneiform-linux-1.0.0/install/lib64  -lcuneiform -lrcorrkegl -lrfrmt -lrmarker -lrblock -lrneg -lrout -lced -lrpic -lrselstr -lrstuff -lrimage -lrline -lrshelllines -lrverline -lcimage -lcfio -lcpage -llns32 -lrdib -lsmetric -lexc -lloc32 -lrreccom -lrpstr -lrstr -lcline -lrcutp -lpass2 -lrbal -lrsadd -lleo32 -levn32 -lfon32 -lctb32 -lmsk32 -ldif32 -lcpu32 -lr3532 -lmmx32 -lrling -lrlings -lcstr -lccom -lstd32 -lcfcompat
endif

CPP=g++ -fopenmp -g -O3
LD=g++ -fopenmp -g -O3 -fPIC
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

ifeq ($(TESSERACT_ENABLE),yes)
CONDITIONAL_DEFINES_TESSERACT=-DTESSERACT_ENABLE
endif

ifeq ($(CUNEIFORM_ENABLE),yes)
CONDITIONAL_DEFINES_CUNEIFORM=-DCUNEIFORM_ENABLE
endif

CONDITIONAL_DEFINES=$(CONDITIONAL_DEFINES_TESSERACT) $(CONDITIONAL_DEFINES_CUNEIFORM)

CPPFLAGS= -g -O2 -fPIC $(OCRADINC) $(MINGWINC) -D_LIB -D_MT $(CONDITIONAL_DEFINES) -Wall $(POTRACEINC) $(GOCRINC) $(MOL_BACKEND_INC) $(TCLAPINC) $(MAGIKINC) $(TESSERACTINC) $(CUNEIFORMINC)

LIBS=$(POTRACELIB) $(OCRADLIB) -lm  $(MAGIKLIB) $(GOCRLIB) $(MOL_BACKEND_LIB)  $(TESSERACTLIB) $(CUNEIFORMLIB)
OBJ = osra.o osra_anisotropic.o osra_ocr.o $(MOL_BACKEND_OBJ) $(MCDLUTIL) unpaper.o


all:	$(OBJ)
	${LD} $(LDFLAGS_STATIC)  -o osra $(OBJ) $(LIBS)


osra.o:	osra.cpp osra.h #pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) $(CPPFLAGS) -c osra.cpp

$(MOL_BACKEND_OBJ): $(MOL_BACKEND_CPP) osra.h
	$(CPP) $(CPPFLAGS) -c $(MOL_BACKEND_CPP)

$(MCDLUTIL):	mcdlutil.cpp mcdlutil.h
	$(CPP) $(CPPFLAGS) -c mcdlutil.cpp

osra_anisotropic.o:	osra_anisotropic.cpp osra.h CImg.h greycstoration.h
	$(CPP) $(CPPFLAGS) -c osra_anisotropic.cpp

osra_ocr.o:	osra_ocr.cpp osra.h #pgm2asc.h output.h list.h unicode.h gocr.h
	$(CPP) $(CPPFLAGS) -c osra_ocr.cpp

clean:	
	-$(RM) -f *.o osra #pgm2asc.h output.h list.h unicode.h gocr.h

#pgm2asc.h: $(GOCRSRC)/pgm2asc.h
#	$(CP) $(GOCRSRC)/pgm2asc.h $(SRCDIR)
#output.h: $(GOCRSRC)/output.h
#	$(CP) $(GOCRSRC)/output.h $(SRCDIR)
#gocr.h: $(GOCRSRC)/gocr.h
#	$(CP) $(GOCRSRC)/gocr.h $(SRCDIR)	
#unicode.h: $(GOCRSRC)/unicode.h
#	$(SED) '/INFINITY/d' $(GOCRSRC)/unicode.h >$(SRCDIR)/unicode.h
#list.h: $(GOCRSRC)/list.h
#	$(SED) 's/struct\ list/struct\ list\_s/' $(GOCRSRC)/list.h >$(SRCDIR)/list.h


unpaper.o:	unpaper.cpp unpaper.h
	$(CPP) $(CPPFLAGS) -c unpaper.cpp