dnl
dnl Autoconf configuration for OSRA.
dnl
dnl Use "autoconf" to generate a configure script.
dnl
dnl Follow this link for complete autoconf macro list:
dnl http://www.gnu.org/software/autoconf/manual/html_node/Autoconf-Macro-Index.html
dnl

AC_INIT([osra], [1.3.6])
dnl AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_BUILD
AC_LANG([C++])
AC_CXX_HAVE_STL

if test "${ac_cv_cxx_have_stl}" != "yes"; then
	AC_MSG_ERROR([STL was not found; make sure you have installed libstdc++-dev])
fi

dnl
dnl Checks for programs
dnl

AC_PROG_CXX
AC_PATH_PROG(RM, rm)
AC_PROG_LN_S
AC_PROG_INSTALL

AC_ARG_ENABLE(
	[docs],
	[AC_HELP_STRING([--enable-docs], [Enable generating of documentation (man pages, API reference, etc)])],
	[
		dnl NaturalDocs executable is present as "/usr/bin/naturaldocs" in Debian and "/usr/bin/NaturalDocs" in Fedora: 
		AC_PATH_PROG(NATURALDOCS, naturaldocs)
		if test -z "${NATURALDOCS}"
		then
			AC_PATH_PROG(NATURALDOCS, NaturalDocs)
			if test -z "${NATURALDOCS}"
			then
				AC_MSG_ERROR([NaturalDocs executable was not found; make sure you have installed NaturalDocs package])
			fi
		fi
		
		AC_PATH_PROG(XSLTPROC, xsltproc)
		if test -z "${XSLTPROC}"
		then
			AC_MSG_ERROR([xsltproc executable was not found; make sure you have installed xsltproc package])
		fi
	]
)

dnl
dnl Include "/usr/local" tree into the list of default locations
dnl

if test -d "/usr/local/include"; then
	CPPFLAGS="-I/usr/local/include ${CPPFLAGS}"
fi
if test -d "/usr/local/lib"; then
	LDFLAGS="-L/usr/local/lib ${LDFLAGS}"
fi

dnl
dnl Compilation/linking configuration options
dnl

AC_ARG_ENABLE(
	[openmp],
	[AC_HELP_STRING([--enable-openmp], [Enable OpenMP extensions for compilers that support it])],
	[
		AC_OPENMP()
		CPPFLAGS="${CPPFLAGS} ${OPENMP_CFLAGS}";
		CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}";
	]
)

dnl Note: "-Wl,-static" affects only libraries ("-l"), following this option.

AC_ARG_ENABLE(
	[static-linking],
	[AC_HELP_STRING([--enable-static-linking], [Link all libraries statically])],
	[LDFLAGS="${LDFLAGS} -static-libgcc -Wl,-static";]
)

AC_ARG_ENABLE(
	[profiling],
	[AC_HELP_STRING([--enable-profiling], [Include verbose profiling/debugging information into the binary])],
	[LDFLAGS="${LDFLAGS} -pg"; CXXFLAGS="${CXXFLAGS} -Wall -g3 -pg";],
	[INSTALL_PROGRAM="${INSTALL_PROGRAM} -s";]
)

dnl
dnl Define includes paths and libraries paths
dnl

dnl tclap library

AX_PROBE_LIBRARY(
	[tclap],
	[CmdLine.h],
	[auto /usr/local/include/tclap],
	[custom location of CmdLine.h from TCLAP (Templatized C++ Command Line Parser Library), see http://tclap.sourceforge.net/]
)

dnl potrace library

AX_PROBE_LIBRARY(
	[potrace],
	[potracelib.h],
	[auto],
	[custom location of potracelib.h from POTRACE (vector tracing library), see http://potrace.sourceforge.net/]
)

dnl gocr library

AX_PROBE_LIBRARY(
	[gocr],
	[pgm2asc.h],
	[auto /usr/local/include/gocr],
	[custom location of pgm2asc.h from GOCR (optical character recognition library), see http://jocr.sourceforge.net/]
)

dnl ocrad library

AX_PROBE_LIBRARY(
	[ocrad],
	[ocradlib.h],
	[auto /usr/local/include],
	[custom location of ocradlib.h from OCRAD (optical character recognition program), see http://www.gnu.org/software/ocrad/ocrad.html]
)

dnl openbabel library

AX_PROBE_LIBRARY(
	[openbabel],
	[openbabel/mol.h openbabel/obconversion.h openbabel/builder.h],
	[/usr/include/openbabel-2.0 /usr/local/include/openbabel-2.0],
	[custom location of openbabel/builder.h from OpenBabel (open source chemistry toolbox), see http://openbabel.sourceforge.net/wiki/Main_Page]
)

dnl
dnl Checks for header files, libraries, types and functions
dnl

dnl limits.h header

osra_header_limits=yes

AC_CHECK_HEADER([limits.h], [], [osra_header_limits=no])
AC_TRY_COMPILE([#include <limits.h>], [ int a = INT_MAX; ], [], [osra_header_limits=no])

if test "${osra_header_limits}" != "yes"; then
	AC_MSG_ERROR([limits.h was not found, do you have libc-dev installed?])
fi

dnl float.h header

osra_header_float=yes

AC_CHECK_HEADER([float.h], [], [osra_header_float=no])
AC_TRY_COMPILE([#include <float.h>], [ float f = FLT_MAX; ], [], [osra_header_float=no])

if test "${osra_header_float}" != "yes"; then
	AC_MSG_ERROR([float.h was not found])
fi

dnl libgen.h header

osra_header_libgen=yes

AC_CHECK_HEADER([libgen.h], [], [osra_header_libgen=no])
AC_TRY_COMPILE([#include <libgen.h>], [ char *a = dirname(""); ], [], [osra_header_libgen=no])

if test "${osra_header_libgen}" != "yes"; then
	AC_MSG_ERROR([libgen.h was not found, do you have libc-dev installed?])
fi

dnl math library

osra_lib_math=yes

AC_CHECK_HEADER([math.h], [], [osra_lib_math=no])
AX_TRY_LINK([m], [#include <math.h>], [ sqrt(1); fabs(1); ], [], [osra_lib_math=no])

if test "${osra_lib_math}" != "yes"; then
	AC_MSG_ERROR([math.h or libm was not found, do you have libc-dev installed?])
fi

dnl tclap library

osra_lib_tclap=yes

AC_TRY_COMPILE([#include <tclap/CmdLine.h>], [
	TCLAP::CmdLine cmd("");
	TCLAP::ValueArg<int> a1("", "", "", false, 0, "");
	cmd.add(a1);
	TCLAP::SwitchArg a2("", "", "", false);
	cmd.add(a2);
	cmd.parse(NULL, NULL);
], [], [osra_lib_tclap=no])

if test "${osra_lib_tclap}" != "yes"; then
	AC_MSG_ERROR([CmdLine.h was not found (use --with-tclap to define custom header location or check config.log)])
fi

dnl potrace library

osra_lib_potrace=yes

AX_TRY_LINK([potrace], [
	extern "C" {
	#include <potracelib.h>
	}
], [
	potrace_state_t state;
	potrace_param_t param;
	potrace_bitmap_t bitmap;
	potrace_path_t path;
	potrace_dpoint_t dpoint;

	potrace_trace(&param, &bitmap);
	potrace_state_free(&state);
	potrace_param_default();
], [], [osra_lib_potrace=no])

if test "${osra_lib_potrace}" != "yes"; then
	AC_MSG_ERROR([potracelib.h or libpotrace was not found (use --with-potrace to define custom header location or check config.log)])
fi

dnl netpbm library check (which is used by GOCR)
dnl (copied from GOCR configure.in)

AC_SEARCH_LIBS(pnm_readpnminit, [netpbm pnm])

dnl gocr library

osra_lib_gocr=yes

dnl I have to add missing global JOB declaration because of the following:
dnl (box.o): In function `cut_box':
dnl (box.o): In function `merge_boxes':
dnl (pixel.o): In function `pixel_filter_by_tree':
dnl (pixel.o): In function `getpixel':
dnl undefined reference to `JOB'

AX_TRY_LINK([Pgm2asc], [
	extern "C" {
	#include <gocr/pgm2asc.h>
	#include <gocr/gocr.h>
	}
	job_t *JOB;
], [
	job_t job;

	job_init(&job);
	job_free(&job);
	pgm2asc(&job);
], [], [osra_lib_gocr=no])

if test "${osra_lib_gocr}" != "yes"; then
	AC_MSG_ERROR([pgm2asc.h or pgm2asc.a was not found (use --with-gocr to define custom header location or check config.log)])
fi

dnl ocrad library

osra_lib_ocrad=yes

AX_TRY_LINK([ocrad], [
	#include <ocradlib.h>
	extern "C" {
	#include <gocr/pgm2asc.h>
	#include <gocr/gocr.h>
	}
	job_t *JOB;
], [
	struct OCRAD_Pixmap opix;
	
	OCRAD_Descriptor * const ocrdes = OCRAD_open();

	OCRAD_get_errno(ocrdes);
	OCRAD_set_image(ocrdes, &opix, 0);
	OCRAD_scale(ocrdes, 2);
	OCRAD_recognize(ocrdes, 0);
	OCRAD_result_first_character(ocrdes);
	OCRAD_close(ocrdes);
], [], [osra_lib_ocrad=no])

if test "${osra_lib_ocrad}" != "yes"; then
	AC_MSG_ERROR([ocradlib.h or libocrad.a was not found (use --with-ocrad to define custom header location)])
fi

dnl openbabel library

osra_lib_openbabel=yes

AX_TRY_LINK([openbabel], [
	#include <openbabel/mol.h>
	#include <openbabel/obconversion.h>
	#include <openbabel/builder.h>
	extern "C" {
	#include <gocr/pgm2asc.h>
	#include <gocr/gocr.h>
	}
	job_t *JOB;
], [
	OpenBabel::OBPairData label;
	OpenBabel::OBAtom atom;
	OpenBabel::OBMol mol;
	OpenBabel::OBBond bond;
	OpenBabel::OBRing ring;
	OpenBabel::OBConversion conversion;

	label.SetAttribute("");
	label.SetValue("");

	atom.IsHeteroatom();
	atom.IsCarbon();
	atom.IsNitrogen();
	atom.IsOxygen();
	atom.IsSulfur();
	atom.IsChiral();
	atom.GetAtomicNum();
	atom.GetFormalCharge();
	atom.SetVector(0, 0, 0);
	atom.SetAtomicNum(0);
	atom.SetFormalCharge(0);
	atom.SetIdx(0);
	atom.Clear();

	mol.NumAtoms();
	mol.NumBonds();
	mol.NumRotors();
	mol.GetAtom(0);
	mol.GetBond(0);
	mol.GetFirstAtom();
	mol.GetSSSR();
	mol.AddAtom(atom);
	mol.AddBond(0, 0, 0, 0);
	mol.CreateAtom();
	mol.AssignSpinMultiplicity(false);
	mol.SetDimension(0);
	mol.SetData(&label);
	std::vector< std::vector<int> > v;
	mol.ContigFragList(v);
	mol.ConnectTheDots();
	mol.Separate();
	mol.StripSalts(0);
	mol.BeginModify();
	mol.EndModify();
	mol.FindRingAtomsAndBonds();

	bond.GetBO();
	bond.GetFlags();
	bond.SetBondOrder(0);
	bond.IsInRing();
	bond.IsDouble();
	bond.IsTriple();
	bond.UnsetUp();
	bond.UnsetDown();
	bond.UnsetAromatic();
	bond.UnsetHash();
	bond.UnsetWedge();
	bond.GetBeginAtom();

	ring.IsAromatic();

	conversion.SetInFormat("");
	conversion.SetOutFormat("");
	conversion.ReadString(&mol, "");
	conversion.Read(&mol);
	conversion.WriteString(&mol, false);

], [], [osra_lib_openbabel=no])

if test "${osra_lib_openbabel}" != "yes"; then
	AC_MSG_ERROR([openbabel headers or library was not found (use --with-openbabel to define custom header location)])
fi

dnl GraphicsMagick library

AC_PATH_PROG(GRAPHICS_MAGICK_CONFIG, GraphicsMagick++-config)

if test "${GRAPHICS_MAGICK_CONFIG}" == ""; then
	AC_MSG_ERROR([GraphicsMagick++-config was not found; make sure you have installed libgraphicsmagick++-dev])
fi

CPPFLAGS="${CPPFLAGS} `${GRAPHICS_MAGICK_CONFIG} --cppflags` -Dcimg_display_type=0"
LDFLAGS="${LDFLAGS} `${GRAPHICS_MAGICK_CONFIG} --ldflags`"
# | sed -r 's#-L/usr/lib(\s+|$)##g'`"

LIBS="${LIBS} `${GRAPHICS_MAGICK_CONFIG} --libs`" 

dnl On some reason these libraries are not included into linker search path as transient dependencies of libGraphicsMagick++ under CygWin. 
if test $build_os == "cygwin"; then
	LIBS="${LIBS} -lX11"
fi

osra_lib_graphicsmagick=yes

AC_CHECK_HEADERS([magick/api.h Magick++.h], [], [osra_lib_graphicsmagick=no])

AX_TRY_LINK([GraphicsMagick++], [
	#include <Magick++.h>
	extern "C" {
	#include <gocr/pgm2asc.h>
	#include <gocr/gocr.h>
	}
	job_t *JOB;
	], [
	Magick::Image image;
	Magick::Color color;
	Magick::ColorRGB colorRGB;
	Magick::ColorGray colorGray;

	image.getPixels(0, 0, 0, 0);
	image.pixelColor(0, 0);
	image.columns();
	image.rows();
	image.matte();
	image.backgroundColor(color);
	image.rotate(0);

	color.alpha();

	colorRGB.red();
	colorRGB.green();
	colorRGB.blue();

	colorGray.shade();
], [], [osra_lib_graphicsmagick=no])

if test "${osra_lib_graphicsmagick}" != "yes"; then
	AC_MSG_ERROR([magick/api.h was not found; make sure you have installed libgraphicsmagick++-dev and all dependent -dev packages])
fi

dnl
dnl Generate Makefile and config.h
dnl

dnl Substitute config.h paths correctly:
test "$prefix" = "NONE" && prefix=$ac_default_prefix

dnl This will resolve 2st-level reference (${datarootdir} and ${PACKAGE_TARNAME}):
datadir=`eval echo "${datadir}"`
docdir=`eval echo "${docdir}"`
dnl This will resolve 1st-level reference (${prefix}):
resolved_datadir=`eval echo "${datadir}"`

dnl C-code needs everything completely resolved:
AC_DEFINE_UNQUOTED(DATA_DIR, ["$resolved_datadir"])
AC_SUBST([resolved_datadir], [$resolved_datadir])

AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_FILES([Makefile.inc doc/manual.sgml])
AC_OUTPUT
AC_MSG_NOTICE([Now you can run "make all install" to compile and install the application.])
