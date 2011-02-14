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
endif

ifdef OSRA_ANDROID
OBJ_ANDROID	:= $(OBJ_LIB) osra_cli.o
endif

all:
	$(RM) -f $(OBJ_CLI)
	$(MAKE) osra
	$(RM) -f $(OBJ_LIB)
	$(MAKE) libosra.a
	$(RM) -f $(OBJ_LIB)
	$(MAKE) libosra.so
ifdef OSRA_JAVA
	$(RM) -f $(OBJ_JAVA)
	$(MAKE) libosra_java.so
endif
ifdef OSRA_ANDROID
	$(RM) -f $(OBJ_ANDROID)
	$(MAKE) libosra_andriod.so
endif
#	Not confusing install rule:
	touch osra libosra.a libosra.so

osra: $(OBJ_CLI)
	$(LD) $(LDFLAGS) -o $@ $(OBJ_CLI) $(LIBS)

libosra.a: CXXFLAGS += -DOSRA_LIB
libosra.a: $(OBJ_LIB)
	$(AR) cru $@ $(OBJ_LIB)
	$(RANLIB) $@

libosra.so: CXXFLAGS += -fPIC -DOSRA_LIB
libosra.so: $(OBJ_LIB)
	$(LD) $(LDFLAGS) -fPIC -shared -Wl,-soname,$@ -o $@ $(OBJ_LIB) $(LIBS)

ifdef OSRA_JAVA
libosra_java.so: CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_JAVA
libosra_java.so: $(OBJ_JAVA)
	$(LD) $(LDFLAGS) -fPIC -shared -Wl,-soname,$@ -o $@ $(OBJ_JAVA) $(LIBS)
endif

ifdef OSRA_ANDROID
libosra_andriod.so: CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_ANDROID
libosra_andriod.so: $(OBJ_ANDROID)
	$(LD) $(LDFLAGS) -fPIC -shared -Wl,-soname,$@ -o $@ $(OBJ_ANDROID) $(LIBS)
endif

Makefile.dep: $(wildcard *.cpp)
	$(CXX) $(CPPFLAGS) -MM $^ > Makefile.dep 

install: osra libosra.a libosra.so
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(libdir) $(DESTDIR)$(includedir)
	$(INSTALL_PROGRAM) osra $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) libosra.so $(DESTDIR)$(libdir)/libosra.so.$(LIB_API_VERSION)
	$(INSTALL_PROGRAM) libosra.a $(DESTDIR)$(libdir)
	$(INSTALL_DATA) osra_lib.h $(DESTDIR)$(includedir)

uninstall:
	$(RM) -f \
		$(DESTDIR)$(bindir)/osra \
		$(DESTDIR)$(libdir)/libosra.so.$(LIB_API_VERSION) \
		$(DESTDIR)$(libdir)/libosra.a \
		$(DESTDIR)$(includedir)/osra_lib.h

clean:
	$(RM) -f *.o osra libosra*.*

distclean: clean
	$(RM) -f config.h Makefile.dep

../Makefile.inc: ../Makefile.inc.in ../config.status
	cd .. && ./config.status
