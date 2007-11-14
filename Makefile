######################### GCC Setup ###################
RES = osra
X11LIBS=-L/usr/X11R6/lib
#X11LIBS=-L/usr/X11R6/lib64
IMLIBS= $(X11LIBS) -lMagick++ -lWand -lMagick -ltiff -lfreetype -ljpeg -lgs -lXext -lSM -lICE -lX11 -lXt -lbz2 -lz
POTRACELIB=-L../../potrace-1.7/src/ -lpotrace
POTRACEINC=-I../../potrace-1.7/src/
GOCRINC= -I../../gocr-0.43/src/ -I../../gocr-0.43/include/
GOCRLIB= -L../../gocr-0.43/src/ -lPgm2asc -lnetpbm
#OPENBABELLIB=-lopenbabel
#OPENBABELINC=-I/usr/include/openbabel-2.0/
OPENBABELLIB=-L/usr/local/lib/ -lopenbabel
OPENBABELINC=-I/usr/local/include/openbabel-2.0/
TCLAPINC=-I/usr/local/include/tclap/
EVG_THIN=evg-thin
OCRADOBJ = arg_parser.o common.o rational.o rectangle.o track.o ucs.o \
       page_image.o page_image_io.o page_image_layout.o \
       bitmap.o block.o profile.o feats.o feats_test0.o feats_test1.o \
       character.o character_r11.o character_r12.o character_r13.o \
       textline.o textline_r2.o textblock.o textpage.o main.o

CPP = g++ -O2  -Wall -DHAVE_CONFIG_H $(POTRACEINC) $(GOCRINC) $(OPENBABELINC) $(TCLAPINC)
LD=gcc
LIBS=$(POTRACELIB) -lm  $(IMLIBS) $(GOCRLIB) $(OPENBABELLIB) 
OBJ = osra.o osra_mol.o $(EVG_THIN).o $(OCRADOBJ)


all:	$(OBJ)
	${LD}  -o ${RES} $(OPTS) $(OBJ) $(LIBS)


osra.o:	osra.cpp osra.h
	$(CPP) -c osra.cpp

osra_mol.o: osra_mol.cpp osra.h
	    $(CPP) -c osra_mol.cpp

${EVG_THIN}.o:	${EVG_THIN}.cc
	$(CPP) -c ${EVG_THIN}.cc

clean:	
	rm -f *.o ${RES}


%.o : %.cc
	$(CPP)  -c -o $@ $<

$(OCRADOBJ)             : Makefile bitmap.h block.h common.h rational.h rectangle.h ucs.h
arg_parser.o        : arg_parser.h
character.o         : character.h
character_r11.o     : character.h profile.h feats.h
character_r12.o     : character.h profile.h feats.h
character_r13.o     : character.h profile.h feats.h
feats.o             : profile.h feats.h
feats_test0.o       : profile.h feats.h
feats_test1.o       : profile.h feats.h
main.o              : arg_parser.h page_image.h textpage.h
page_image.o        : page_image.h
page_image_io.o     : page_image.h
page_image_layout.o : track.h page_image.h
profile.o           : profile.h
textblock.o         : track.h character.h page_image.h textline.h textblock.h
textline.o          : track.h character.h page_image.h textline.h
textline_r2.o       : track.h character.h textline.h
textpage.o          : track.h character.h page_image.h textline.h textblock.h textpage.h
track.o             : track.h

