#
# This makefile is responsible for building the executable.
#

OBJ		= osra.o osra_anisotropic.o osra_ocr.o osra_openbabel.o mcdlutil.o unpaper.o

include ../Makefile.inc
include Makefile.dep

all: osra

osra: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

Makefile.dep: $(patsubst %.o,%.cpp,$(OBJ))
	$(CXX) $(CPPFLAGS) -MM $? > Makefile.dep 

install: osra
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) osra $(DESTDIR)$(bindir)

uninstall:
	$(RM) -f $(DESTDIR)$(bindir)/osra

clean:
	$(RM) -f *.o osra

distclean: clean
	$(RM) -f config.h Makefile.dep

#../Makefile.inc: ../Makefile.inc.in ../config.status
#	cd .. && ./config.status
