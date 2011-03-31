#
# This makefile is responsible for building the executable.
#

include ../Makefile.inc
include Makefile.dep

.PHONY: clean_obj

LIB_VERSION	:= $(LIB_MAJOR_VERSION).$(LIB_MINOR_VERSION).$(LIB_PATCH_VERSION)

TARGETS		:= osra

OBJ_LIB		:= osra.o osra_anisotropic.o osra_ocr.o osra_openbabel.o mcdlutil.o unpaper.o 

ifdef TESSERACT_LIB
OBJ_LIB		+= osra_ocr_tesseract.o
endif

OBJ_CLI		:= $(OBJ_LIB) osra_cli.o

ifdef OSRA_LIB
TARGETS		+= libosra.a libosra.so
endif

ifdef OSRA_JAVA
OBJ_JAVA	:= $(OBJ_LIB) osra_java.o
TARGETS		+= libosra_java.so
endif

ifdef OSRA_ANDROID
OBJ_ANDROID	:= $(OBJ_LIB) osra_cli.o
endif

all:
	$(RM) -f $(OBJ_CLI)
	$(MAKE) osra
ifdef OSRA_LIB
	$(RM) -f $(OBJ_LIB)
	$(MAKE) libosra.a
	$(RM) -f $(OBJ_LIB)
	$(MAKE) libosra.so
endif
ifdef OSRA_JAVA
	$(RM) -f $(OBJ_JAVA)
	$(MAKE) libosra_java.so
endif
ifdef OSRA_ANDROID
	$(RM) -f $(OBJ_ANDROID)
	$(MAKE) libosra_andriod.so
endif
#	We have to update the timestaps of the targets, otherwise "install" target will try to re-link and will cause missed symbols:
	touch $(TARGETS)

osra: $(OBJ_CLI)
	$(LD) $(LDFLAGS) -o $@ $(OBJ_CLI) $(LIBS)

ifdef OSRA_LIB
libosra.a: CXXFLAGS += -DOSRA_LIB
libosra.a: $(OBJ_LIB)
	$(AR) cru $@ $(OBJ_LIB)
	$(RANLIB) $@

libosra.so: CXXFLAGS += -fPIC -DOSRA_LIB
libosra.so: $(OBJ_LIB)
	$(LD) $(LDFLAGS) -fPIC -shared -Wl,-soname,$@.$(LIB_MAJOR_VERSION) -o $@ $(OBJ_LIB) $(LIBS)
endif

ifdef OSRA_JAVA
libosra_java.so: CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_JAVA
libosra_java.so: $(OBJ_JAVA)
	$(LD) $(LDFLAGS) -fPIC -shared -Wl,-soname,$@.$(LIB_MAJOR_VERSION) -o $@ $(OBJ_JAVA) $(LIBS)
endif

ifdef OSRA_ANDROID
libosra_andriod.so: CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_ANDROID
libosra_andriod.so: $(OBJ_ANDROID)
	$(LD) $(LDFLAGS) -fPIC -shared -Wl,-soname,$@ -o $@ $(OBJ_ANDROID) $(LIBS)
endif

Makefile.dep: $(patsubst %.o,%.cpp,$(OBJ_CLI))
	$(CXX) $(CPPFLAGS) -MM $^ > Makefile.dep 

install: $(TARGETS)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) osra $(DESTDIR)$(bindir)
ifdef OSRA_LIB
	$(INSTALL_DIR) $(DESTDIR)$(libdir) $(DESTDIR)$(includedir) $(DESTDIR)$(libdir)/pkgconfig
	$(INSTALL_PROGRAM) libosra.so $(DESTDIR)$(libdir)/libosra.so.$(LIB_VERSION)
	$(LN_S) libosra.so.$(LIB_VERSION) $(DESTDIR)$(libdir)/libosra.so.$(LIB_MAJOR_VERSION)
	$(INSTALL_PROGRAM) libosra.a $(DESTDIR)$(libdir)
	$(INSTALL_DATA) osra_lib.h $(DESTDIR)$(includedir)
	$(INSTALL_DATA) ../package/linux/osra.pc $(DESTDIR)$(libdir)/pkgconfig
endif
ifdef OSRA_JAVA
	$(INSTALL_DIR) $(DESTDIR)$(libdir)
	$(INSTALL_PROGRAM) libosra_java.so $(DESTDIR)$(libdir)/libosra_java.so.$(LIB_VERSION)
	$(LN_S) libosra_java.so.$(LIB_VERSION) $(DESTDIR)$(libdir)/libosra_java.so.$(LIB_MAJOR_VERSION)
endif

# "install" and "make" tools auomatically autodetect the extension for executables, but "rm" needs extension correction:
uninstall:
	$(RM) -f \
		$(DESTDIR)$(bindir)/osra$(EXEEXT) \
		$(DESTDIR)$(libdir)/libosra.so.$(LIB_VERSION) \
		$(DESTDIR)$(libdir)/libosra.a \
		$(DESTDIR)$(libdir)/libosra_java.so.$(LIB_VERSION) \
		$(DESTDIR)$(includedir)/osra_lib.h \
		$(DESTDIR)$(libdir)/pkgconfig/osra.pc

clean:
	$(RM) -f *.o osra$(EXEEXT) libosra*.*

distclean: clean
	$(RM) -f config.h Makefile.dep

../Makefile.inc: ../Makefile.inc.in ../config.status
	cd .. && ./config.status
