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
