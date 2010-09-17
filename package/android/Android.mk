LOCAL_PATH := $(call my-dir)

include $(CLEAR_VARS)

LOCAL_MODULE := osra

LOCAL_CFLAGS := -DANDROID

LOCAL_SRC_FILES := osra.cpp osra_ocr.cpp osra_openbabel.cpp osra_anisotropic.cpp mcdlutil.cpp unpaper.cpp

LOCAL_C_INCLUDES := /home/igor/ndk-crystax/potrace-1.8/jni/
LOCAL_C_INCLUDES += /home/igor/ndk-crystax/openbabel-2.2.3/include
LOCAL_C_INCLUDES += /home/igor/ndk-crystax/GraphicsMagick-1.3.8/Magick++/jni/ /home/igor/ndk-crystax/GraphicsMagick-1.3.8/
LOCAL_C_INCLUDES += /home/igor/ndk-crystax/tclap-1.1.0/include/tclap/ /home/igor/ndk-crystax/tclap-1.1.0/include
LOCAL_C_INCLUDES += /home/igor/ndk-crystax/gocr-0.48-patched/jni/
LOCAL_C_INCLUDES += /home/igor/ndk-crystax/ocrad-0.21-pre1/jni/



LOCAL_LDLIBS := -L/home/igor/ndk-crystax/potrace-1.8/obj/local/armeabi/ -lpotrace
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/openbabel-2.2.3/libs/armeabi/ -lopenbabel
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/GraphicsMagick-1.3.8/Magick++/obj/local/armeabi/ -lGraphicsMagick++
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/GraphicsMagick-1.3.8/magick/obj/local/armeabi/ -lGraphicsMagick
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/gocr-0.48-patched/obj/local/armeabi/ -lPgm2acs
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/ocrad-0.21-pre1/obj/local/armeabi -locrad
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/libjasper/obj/local/armeabi/ -ljasper
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/jpeg-8b/obj/local/armeabi/ -ljpeg
LOCAL_LDLIBS += -L/home/igor/ndk-crystax/libpng-1.4.3/obj/local/armeabi/ -lpng
LOCAL_LDLIBS += -lm -lz

include $(BUILD_SHARED_LIBRARY)

