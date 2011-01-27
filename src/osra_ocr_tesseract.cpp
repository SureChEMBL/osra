#include <cstdio>
#include <iostream>
#include <cstdlib>
#include "config.h"
#ifdef HAVE_TESSERACT_LIB

#include <tesseract/baseapi.h>

char *tesseract_backend(unsigned char *p, int x1, int y1, int x2, int y2)
{
  tesseract::TessBaseAPI tess;
  tess.Init(NULL, "eng", NULL, 0, false);
  char* text = tess.TesseractRect(p, 1, x2 - x1 + 1, 0, 0, x2 - x1 + 1, y2 - y1 + 1);
  tess.End();

  return(text);
}

#endif
