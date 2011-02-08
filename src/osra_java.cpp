/******************************************************************************
 OSRA: Optical Structure Recognition

 This is a U.S. Government work (2007-2010) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2010 (igorf@helix.nih.gov)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

#include <jni.h>
#include <stdlib.h> // calloc(), free()

#include <string> // std::string
#include <ostream> // std:ostream
#include <sstream> // std:stringstream

#include "config.h" // PACKAGE_VERSION

using namespace std;

#ifdef ANDROID
int main(int argc, char **argv, const char *image_data, int image_length, ostream &output_structure_stream);

extern "C" {
  JNIEXPORT jstring JNICALL Java_cadd_osra_main_runosra_nativeosra(JNIEnv *j_env, jobject j_this, jobjectArray j_arr,
      jbyteArray j_image_data);
};

JNIEXPORT jstring JNICALL Java_cadd_osra_main_runosra_nativeosra(JNIEnv *j_env, jobject j_this, jobjectArray j_arr,
    jbyteArray j_image_data)
{
  int argc = j_env->GetArrayLength(j_arr);
  char **argv = (char **) calloc(argc, sizeof(char*));

  // Convert from Java UTF-16 string to UTF-8 string (assuming that all chars are ASCII):
  for (int i = 0; i < argc; i++)
    {
      jstring j_str = (jstring) j_env->GetObjectArrayElement(j_arr, i);
      argv[i] = j_env->GetStringUTFChars(j_str, NULL);
    }

  char *image_data = (char *) j_env->GetByteArrayElements(j_image_data, NULL);

  int result = 0;
  stringstream output_structure_stream;

  if (image_data != NULL)
    {
      result = main(argc, argv, image_data, j_env->GetArrayLength(j_image_data), output_structure_stream);

      j_env->ReleaseByteArrayElements(j_image_data, (jbyte *) image_data, JNI_ABORT);
    }

  // Release resources:
  for (int i = 0; i < argc; i++)
    {
      jstring j_str = (jstring) j_env->GetObjectArrayElement(j_arr, i);
      j_env->ReleaseStringUTFChars(j_str, argv[i]);
    }

  free(argv);

  if (result != 0)
    return j_env->NewStringUTF("");

  return j_env->NewStringUTF(output_structure_stream.str().c_str());
}
#endif

#ifdef OSRA_JAVA
#include "osra_lib.h"

extern "C" {
  /*
   * Class:     net_sourceforge_osra_OsraLib
   * Method:    processImage
   * Signature: ([BLjava/io/Writer;Ljava/lang/String;ZZZ)I
   */
  JNIEXPORT jint JNICALL Java_net_sourceforge_osra_OsraLib_processImage(JNIEnv *, jobject, jbyteArray, jobject, jstring, jboolean, jboolean, jboolean);

  /*
   * Class:     net_sourceforge_osra_OsraLib
   * Method:    getVersion
   * Signature: ()Ljava/lang/String;
   */
  JNIEXPORT jstring JNICALL Java_net_sourceforge_osra_OsraLib_getVersion(JNIEnv *, jobject);
}

JNIEXPORT jint JNICALL Java_net_sourceforge_osra_OsraLib_processImage(JNIEnv *j_env, jobject j_this, jbyteArray j_image_data, jobject j_writer, jstring j_format,
                                                                         jboolean j_output_confidence, jboolean j_output_coordinates,
                                                                         jboolean j_output_avg_bond_length)
{
  const char *format = j_env->GetStringUTFChars(j_format, NULL);
  const char *image_data = (char *) j_env->GetByteArrayElements(j_image_data, NULL);

  int result = -1;

  if (image_data != NULL)
    {
      stringstream output_structure_stream;

      result = osra_process_image(
                                  image_data,
                                  j_env->GetArrayLength(j_image_data),
                                  output_structure_stream,
                                  0,
                                  false,
                                  0,
                                  0,
                                  0,
                                  false,
                                  string(format),
                                  j_output_confidence,
                                  false,
                                  true,
                                  j_output_coordinates,
                                  j_output_avg_bond_length
      );

      j_env->ReleaseByteArrayElements(j_image_data, (jbyte *) image_data, JNI_ABORT);

      // Locate java.io.Writer#write(char[]) method:
      jclass j_writer_class = j_env->FindClass("java/io/Writer");
      jmethodID write_method_id = j_env->GetMethodID(j_writer_class, "write", "([C)V");

      // Create a char[] that holds the string characters:
      jsize len = output_structure_stream.str().length();
      jcharArray j_char_arr = j_env->NewCharArray(len);
      j_env->SetCharArrayRegion(j_char_arr, 0, len, (jchar *) output_structure_stream.str().c_str());

      j_env->CallVoidMethod(j_writer, write_method_id, j_char_arr);

      j_env->DeleteLocalRef(j_writer_class);
      j_env->DeleteLocalRef(j_char_arr);

      jstring output_structure = j_env->NewStringUTF(output_structure_stream.str().c_str());
    }

  j_env->ReleaseStringUTFChars(j_format, format);

  return result;
}

JNIEXPORT jstring JNICALL Java_net_sourceforge_osra_OsraLib_getVersion(JNIEnv *j_env, jobject j_this)
{
  return j_env->NewStringUTF(PACKAGE_VERSION);
}
#endif
