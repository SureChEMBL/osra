OSRA_DIR := $(call my-dir)/../
#  potrace: Need to copy config.h to src/
POTRACE_DIR := $(OSRA_DIR)/../potrace-1.8/
# openbabel:  Need to copy babelconfig.h to include/openbabel; Copy InChI headers: from include/inchi102/*.h to src/formats/inchi102/
OPENBABEL_DIR := $(OSRA_DIR)/../openbabel-2.2.3/
# gocr:  Need to copy config.h to src/
GOCR_DIR := $(OSRA_DIR)/../gocr-0.48-patched/
OCRAD_DIR := $(OSRA_DIR)/../ocrad-0.21-pre1/
TCLAP_DIR := $(OSRA_DIR)/../tclap-1.1.0/
# Magick++: Need to copy magick_config.h to magick/, magick: Edit coders/png.c to comment out lines 95-108 (PNG_LIBPNG_VER < 10400)
GRAPHICSMAGICK_DIR := $(OSRA_DIR)/../GraphicsMagick-1.3.8/
LIBPNG_DIR := $(OSRA_DIR)/../libpng-1.4.3/
# jpeg: Need to copy jconfig.h to top level
JPEG_DIR := $(OSRA_DIR)/../jpeg-8b/
OSRA_SRC_DIR := $(OSRA_DIR)/src/

LOCAL_PATH := $(OSRA_SRC_DIR)
include $(CLEAR_VARS)

LOCAL_MODULE := osra

LOCAL_CFLAGS := -DANDROID

LOCAL_SRC_FILES := osra.cpp osra_ocr.cpp osra_openbabel.cpp osra_anisotropic.cpp mcdlutil.cpp unpaper.cpp

LOCAL_C_INCLUDES := $(POTRACE_DIR)/src/
LOCAL_C_INCLUDES += $(OPENBABEL_DIR)/include
LOCAL_C_INCLUDES += $(GRAPHICSMAGICK_DIR)/Magick++/lib/ $(GRAPHICSMAGICK_DIR)
LOCAL_C_INCLUDES += $(TCLAP_DIR)/include/tclap/ $(TCLAP_DIR)/include/
LOCAL_C_INCLUDES += $(GOCR_DIR)/src/
LOCAL_C_INCLUDES += $(OCRAD_DIR)


LOCAL_LDLIBS := -lz

LOCAL_STATIC_LIBRARIES := potrace GraphicsMagick++ GraphicsMagick Pgm2acs ocrad jpeg png
LOCAL_SHARED_LIBRARIES := openbabel

include $(BUILD_SHARED_LIBRARY)

# potrace: Need to copy config.h to src/
LOCAL_PATH := $(POTRACE_DIR)/src/

include $(CLEAR_VARS)

LOCAL_MODULE    := potrace
LOCAL_SRC_FILES := curve.c trace.c decompose.c potracelib.c
LOCAL_CFLAGS := -DHAVE_CONFIG_H
include $(BUILD_STATIC_LIBRARY)

# openbabel: Need to copy babelconfig.h to include/openbabel; Copy InChI headers: from include/inchi102/*.h to src/formats/inchi102/
LOCAL_PATH := $(OPENBABEL_DIR)/src/
OB_SRC_TOP_DIR := $(LOCAL_PATH)

include $(CLEAR_VARS)

LOCAL_MODULE := openbabel
LOCAL_CFLAGS := -DANDROID

LOCAL_SRC_FILES := atom.cpp alias.cpp base.cpp bitvec.cpp bond.cpp bondtyper.cpp canon.cpp builder.cpp chains.cpp chiral.cpp data.cpp descriptor.cpp fingerprint.cpp forcefield.cpp format.cpp generic.cpp griddata.cpp grid.cpp locale.cpp kekulize.cpp matrix.cpp mol.cpp molchrg.cpp obconversion.cpp oberror.cpp obiter.cpp obmolecformat.cpp obutil.cpp op.cpp ops.cpp parsmart.cpp patty.cpp phmodel.cpp plugin.cpp pointgroup.cpp rand.cpp residue.cpp ring.cpp rotor.cpp rotamer.cpp tokenst.cpp typer.cpp transform.cpp 
LOCAL_SRC_FILES += formats/smilesformat.cpp formats/getinchi.cpp formats/inchiformat.cpp

LOCAL_C_INCLUDES := $(LOCAL_PATH)/../include $(LOCAL_PATH)/../data

LOCAL_STATIC_LIBRARIES := math fingerprints inchi forcefields ops descriptors

include $(BUILD_SHARED_LIBRARY)

# math library
LOCAL_PATH := $(OB_SRC_TOP_DIR)/math/

include $(CLEAR_VARS)

LOCAL_MODULE := math

LOCAL_SRC_FILES := matrix3x3.cpp spacegroup.cpp transform3d.cpp vector3.cpp
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../../include


include $(BUILD_STATIC_LIBRARY)

# fingerprints library
LOCAL_PATH := $(OB_SRC_TOP_DIR)/fingerprints/

include $(CLEAR_VARS)

LOCAL_MODULE := fingerprints

LOCAL_SRC_FILES := finger2.cpp finger3.cpp
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../../include

include $(BUILD_STATIC_LIBRARY)

# inchi library
LOCAL_PATH := $(OB_SRC_TOP_DIR)/formats/inchi102/

include $(CLEAR_VARS)

LOCAL_MODULE := inchi

LOCAL_SRC_FILES :=  ichi_bns.c ichi_io.c ichican2.c ichicano.c ichicans.c ichiisot.c ichilnct.c ichimak2.c ichimake.c ichimap1.c ichimap2.c ichimap4.c ichinorm.c ichiparm.c ichiprt1.c ichiprt2.c ichiprt3.c ichiqueu.c ichiread.c ichiring.c ichirvr1.c ichirvr2.c ichirvr3.c ichirvr4.c ichirvr5.c ichirvr6.c ichirvr7.c ichisort.c ichister.c ichitaut.c ikey_base26.c ikey_dll.c inchi_dll.c inchi_dll_a.c inchi_dll_a2.c inchi_dll_main.c runichi.c sha2.c strutil.c util.c


include $(BUILD_STATIC_LIBRARY)

# forcefields library
LOCAL_PATH := $(OB_SRC_TOP_DIR)/forcefields/

include $(CLEAR_VARS)

LOCAL_MODULE := forcefields

LOCAL_SRC_FILES := forcefieldghemical.cpp forcefieldmmff94.cpp forcefielduff.cpp
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../../include

include $(BUILD_STATIC_LIBRARY)

# ops library
LOCAL_PATH := $(OB_SRC_TOP_DIR)/ops/

include $(CLEAR_VARS)

LOCAL_MODULE := ops

LOCAL_SRC_FILES := addpolarh.cpp fillUC.cpp gen3d.cpp loader.cpp optransform.cpp
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../../include

include $(BUILD_STATIC_LIBRARY)

#descriptors library
LOCAL_PATH := $(OB_SRC_TOP_DIR)/descriptors/

include $(CLEAR_VARS)

LOCAL_MODULE := descriptors

LOCAL_SRC_FILES := cmpdfilter.cpp filters.cpp groupcontrib.cpp smartsdescriptors.cpp
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../../include


include $(BUILD_STATIC_LIBRARY)

# gocr: Need to copy config.h to src/
LOCAL_PATH :=  $(GOCR_DIR)/src/

include $(CLEAR_VARS)

LOCAL_MODULE    := Pgm2acs
LOCAL_SRC_FILES := pgm2asc.c box.c database.c detect.c barcode.c lines.c list.c ocr0.c ocr0n.c ocr1.c otsu.c output.c pixel.c unicode.c remove.c pnm.c pcx.c progress.c job.c
LOCAL_ALLOW_UNDEFINED_SYMBOLS := true
include $(BUILD_STATIC_LIBRARY)

# ocrad 
LOCAL_PATH := $(OCRAD_DIR)

include $(CLEAR_VARS)

LOCAL_CPP_EXTENSION := .cc
LOCAL_MODULE    := ocrad
LOCAL_SRC_FILES := ocradlib.cc common.cc mask.cc rational.cc rectangle.cc track.cc ucs.cc page_image.cc page_image_io.cc bitmap.cc blob.cc profile.cc feats.cc feats_test0.cc feats_test1.cc character.cc character_r11.cc character_r12.cc character_r13.cc textline.cc textline_r2.cc textblock.cc textpage.cc


include $(BUILD_STATIC_LIBRARY)

# libpng
LOCAL_PATH := $(LIBPNG_DIR)

include $(CLEAR_VARS)

LOCAL_MODULE    := png
LOCAL_SRC_FILES := png.c pngset.c pngget.c pngrutil.c pngtrans.c pngwutil.c pngread.c pngrio.c pngwio.c pngwrite.c pngrtran.c pngwtran.c pngmem.c pngerror.c pngpread.c 

include $(BUILD_STATIC_LIBRARY)


# jpeg: Need to copy jconfig.h to top level
LOCAL_PATH := $(JPEG_DIR)

include $(CLEAR_VARS)

LOCAL_MODULE    := jpeg
LOCAL_SRC_FILES := jaricom.c jcapimin.c jcapistd.c jcarith.c jccoefct.c jccolor.c jcdctmgr.c jchuff.c jcinit.c jcmainct.c jcmarker.c jcmaster.c  jcomapi.c jcparam.c jcprepct.c jcsample.c jctrans.c jdapimin.c jdapistd.c jdarith.c jdatadst.c jdatasrc.c jdcoefct.c jdcolor.c jddctmgr.c jdhuff.c jdinput.c jdmainct.c jdmarker.c jdmaster.c jdmerge.c jdpostct.c jdsample.c jdtrans.c jerror.c jfdctflt.c jfdctfst.c jfdctint.c jidctflt.c jidctfst.c jidctint.c jquant1.c jquant2.c jutils.c jmemmgr.c jmemnobs.c

include $(BUILD_STATIC_LIBRARY)

# Magick++: Need to copy magick_config.h to magick/
LOCAL_PATH := $(GRAPHICSMAGICK_DIR)/Magick++/lib/

include $(CLEAR_VARS)
LOCAL_CFLAGS := -DHAVE_CONFIG_H
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../../
LOCAL_MODULE    := GraphicsMagick++
LOCAL_SRC_FILES := Blob.cpp BlobRef.cpp CoderInfo.cpp Color.cpp Drawable.cpp Exception.cpp Functions.cpp Geometry.cpp Image.cpp ImageRef.cpp Montage.cpp Options.cpp Pixels.cpp STL.cpp Thread.cpp TypeMetric.cpp

include $(BUILD_STATIC_LIBRARY)


# magick: Edit coders/png.c to comment out lines 95-108 (PNG_LIBPNG_VER < 10400)
LOCAL_PATH := $(GRAPHICSMAGICK_DIR)/magick/

include $(CLEAR_VARS)
LOCAL_CFLAGS := -DHAVE_CONFIG_H
LOCAL_C_INCLUDES := $(LOCAL_PATH)/../ $(LOCAL_PATH)/../../jpeg-8b/  $(LOCAL_PATH)/../../libpng-1.4.3/
LOCAL_MODULE    := GraphicsMagick
LOCAL_SRC_FILES := analyze.c annotate.c attribute.c average.c bit_stream.c blob.c cdl.c channel.c compare.c confirm_access.c color.c color_lookup.c colormap.c colorspace.c command.c composite.c compress.c constitute.c decorate.c delegate.c deprecate.c describe.c draw.c effect.c enhance.c enum_strings.c error.c fx.c gem.c gradient.c hclut.c image.c list.c locale.c log.c magic.c magick.c magick_endian.c map.c memory.c module.c monitor.c montage.c omp_data_view.c operator.c paint.c pixel_cache.c pixel_iterator.c plasma.c profile.c quantize.c registry.c random.c render.c resize.c resource.c segment.c semaphore.c shear.c signature.c static.c statistics.c tempfile.c texture.c timer.c transform.c tsd.c type.c unix_port.c utility.c version.c ../coders/art.c ../coders/avi.c ../coders/avs.c ../coders/bmp.c ../coders/cals.c ../coders/caption.c ../coders/cineon.c ../coders/cmyk.c ../coders/cut.c ../coders/dcm.c ../coders/dcraw.c ../coders/dib.c ../coders/dpx.c ../coders/fax.c ../coders/fits.c ../coders/gif.c ../coders/gradient.c ../coders/gray.c ../coders/histogram.c ../coders/hrz.c ../coders/html.c ../coders/icon.c ../coders/identity.c ../coders/label.c ../coders/locale.c ../coders/logo.c ../coders/map.c ../coders/mat.c ../coders/matte.c ../coders/meta.c ../coders/miff.c ../coders/mono.c ../coders/mpc.c ../coders/mpeg.c ../coders/mpr.c ../coders/msl.c ../coders/mtv.c ../coders/mvg.c ../coders/null.c ../coders/otb.c ../coders/palm.c ../coders/pcd.c ../coders/pcl.c ../coders/pcx.c ../coders/pdb.c ../coders/pdf.c ../coders/pict.c ../coders/pix.c ../coders/plasma.c ../coders/pnm.c ../coders/preview.c ../coders/ps.c ../coders/ps2.c ../coders/ps3.c ../coders/psd.c ../coders/pwp.c ../coders/rgb.c ../coders/rla.c ../coders/rle.c ../coders/sct.c ../coders/sfw.c ../coders/sgi.c ../coders/stegano.c ../coders/sun.c ../coders/svg.c ../coders/tga.c ../coders/tile.c ../coders/tim.c ../coders/topol.c ../coders/ttf.c ../coders/txt.c ../coders/uil.c ../coders/url.c ../coders/uyvy.c ../coders/vicar.c ../coders/vid.c ../coders/viff.c ../coders/wbmp.c ../coders/wmf.c ../coders/wpg.c ../coders/xbm.c ../coders/xc.c ../coders/xcf.c ../coders/xpm.c ../coders/yuv.c ../coders/jpeg.c ../coders/jp2.c ../coders/png.c ../filters/analyze.c


include $(BUILD_STATIC_LIBRARY)

