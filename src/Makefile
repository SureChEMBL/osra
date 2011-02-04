#
# This makefile is responsible for building the executable.
#

include ../Makefile.inc
include Makefile.dep

.PHONY: clean_obj

LIB_API_VERSION	:= 1.0.0

OBJ_LIB		:= osra.o osra_anisotropic.o osra_ocr.o osra_openbabel.o mcdlutil.o unpaper.o 

ALL			:= osra clean_obj libosra.so clean_obj

ifdef TESSERACT_LIB
OBJ_LIB		+= osra_ocr_tesseract.o
endif

OBJ_CLI		:= $(OBJ_LIB) osra_cli.o

ifdef OSRA_JAVA
OBJ_JAVA	:= $(OBJ_LIB) osra_java.o
ALL			+= libosra_java.so clean_obj
endif

#all: $(ALL)

all:
	$(MAKE) clean_obj osra
	$(MAKE) clean_obj libosra.so
ifdef OSRA_JAVA
	$(MAKE) clean_obj libosra_java.so
endif

osra: $(OBJ_CLI)
	$(LD) $(LDFLAGS) -o $@ $(OBJ_CLI) $(LIBS)

libosra.so: CXXFLAGS += -fPIC -DOSRA_LIB
libosra.so: $(OBJ_LIB)
	$(LD) $(LDFLAGS) -fPIC -shared -o $@ $(OBJ_LIB) $(LIBS)
	mv $@ $@.$(LIB_API_VERSION)

libosra_java.so: CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_JAVA
libosra_java.so: $(OBJ_JAVA)
	$(LD) $(LDFLAGS) -fPIC -shared -o $@ $(OBJ_LIB) $(LIBS)
	mv $@ $@.$(LIB_API_VERSION)

Makefile.dep: $(patsubst %.o,%.cpp,$(OBJ))
	$(CXX) $(CPPFLAGS) -MM $? > Makefile.dep 

install: osra
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) osra $(DESTDIR)$(bindir)

uninstall:
	$(RM) -f $(DESTDIR)$(bindir)/osra

clean_obj:
	$(RM) -f *.o

clean: clean_obj
	$(RM) -f osra libosra*.so*

distclean: clean
	$(RM) -f config.h Makefile.dep

../Makefile.inc: ../Makefile.inc.in ../config.status
	cd .. && ./config.status
