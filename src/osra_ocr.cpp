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

#include <algorithm>
#include <cstdio>
#include <vector>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>

#include "config.h"

extern "C" {
#include <pgm2asc.h>
}

#include <ocradlib.h>

#ifdef HAVE_CUNEIFORM_LIB
#include <cuneiform.h>
#endif

#include "osra.h"

#ifdef HAVE_TESSERACT_LIB
void osra_tesseract_init();
void osra_tesseract_release();
char osra_tesseract_ocr(unsigned char *pixel_map, int x1, int y1, int x2, int y2, string &char_filter);
#endif

// Global GOCR variable (omg) both for 0.48-0.49 and 0.50 versions:
job_t *JOB;
job_t *OCR_JOB;

/**
 * THRESHOLD is the graylevel binarization threshold.
 * dropx and dropy are the coordinates for the starting point from where the connected component (the image of the character) will be searched for.
 * Very often there is no bounding rectangle that would exclude all extraneous pieces without cutting into the character itself.
 * Those pieces confuse the OCR libraries quite a bit, so it's better to extract connected components (all characters that OSRA needs to resolve
 * are luckily single connected components) and leave the extra bits out.
 */

void osra_ocr_init()
{
#ifdef HAVE_CUNEIFORM_LIB
  int langcode = LANG_ENGLISH;
  Bool dotmatrix = 0;
  Bool fax = 0;
  Bool onecolumn = 1;
  PUMA_Init(0, 0);
  PUMA_SetImportData(PUMA_Word32_Language, &langcode);
  PUMA_SetImportData(PUMA_Bool32_DotMatrix, &dotmatrix);
  PUMA_SetImportData(PUMA_Bool32_Fax100, &fax);
  PUMA_SetImportData(PUMA_Bool32_OneColumn, &onecolumn);
#endif

#ifdef HAVE_TESSERACT_LIB
  osra_tesseract_init();
#endif
}

void osra_ocr_release()
{
#ifdef HAVE_CUNEIFORM_LIB
  PUMA_Done();
#endif
#ifdef HAVE_TESSERACT_LIB
  osra_tesseract_release();
#endif
}

//  Function: osra_gocr_ocr()
//    Make an attempt to OCR the image box with GOCR engine.
//
//  Parameters:
//    job_t - includes pixel map and character filter
//
//  Returns:
//    0 in case the recognition failed or valid alphanumeric character
char osra_gocr_ocr(job_t &gocr_job)
{
  OCR_JOB = &gocr_job;
  JOB = &gocr_job;
  try
    {
      pgm2asc(&gocr_job);
    }
  catch (...)
    {
    }

  char *l = (char *) gocr_job.res.linelist.start.next->data;

  if (l != NULL && strlen(l) == 1 && isalnum(l[0]))
    return l[0];

  return 0;
}

//  Function: osra_ocrad_ocr()
//    Make an attempt to OCR the image box with OCRAD engine.
//
//  Parameters:
//    ocrad_pixmap - includes pixel map and the image mode
//    char_filter - character filter
//
//  Returns:
//    0 in case the recognition failed or valid alphanumeric character
char osra_ocrad_ocr(const OCRAD_Pixmap * const ocrad_pixmap, string &char_filter)
{
  char result = 0;
  string line;

  OCRAD_Descriptor * const ocrad_res = OCRAD_open();

  // If the box height is less than 10px, it should be scaled up a bit, otherwise OCRAD is unable to catch it:
  if (ocrad_res && OCRAD_get_errno(ocrad_res) == OCRAD_ok && OCRAD_set_image(ocrad_res, ocrad_pixmap, 0) == 0
      && (ocrad_pixmap->height >= 10 || OCRAD_scale(ocrad_res, 2) == 0) && OCRAD_recognize(ocrad_res, 0) == 0)
    {
      result = OCRAD_result_first_character(ocrad_res);
      if (OCRAD_result_blocks(ocrad_res) >= 1 && OCRAD_result_lines(ocrad_res, 0) && OCRAD_result_line(
            ocrad_res, 0, 0) != 0)
        line = OCRAD_result_line(ocrad_res, 0, 0);
    }

  OCRAD_close(ocrad_res);

  // TODO: Why line should have 0 or 1 characters? Give examples...
  if (line.length() > 2 || !isalnum(result) || char_filter.find(result, 0) == string::npos)
    return 0;

  return result;
}

//  Function: osra_cuneiform_ocr()
//    Make an attempt to OCR the image box with Cuneiform engine.
//
//  Parameters:
//    cuneiform_img - pixel map
//    verbose - if set, then output intermediate results
//    char_filter - character filter
//
//  Returns:
//    0 in case the recognition failed or valid alphanumeric character
char osra_cuneiform_ocr(Magick::Image &cuneiform_img, bool verbose, string &char_filter)
{
  char str[256];

  Magick::Blob blob;
  cuneiform_img.write(&blob, "DIB");
  size_t data_size = blob.length();
  char *dib = new char[data_size];
  memcpy(dib, blob.data(), data_size);

  PUMA_XOpen(dib, NULL);
  PUMA_XFinalRecognition();
  PUMA_SaveToMemory(NULL, PUMA_TOTEXT, PUMA_CODE_ASCII, str, sizeof(str) - 1);
  PUMA_XClose();

  delete[] dib;

  if (verbose)
    cout << "Cuneiform: " << str[0] << "|" << str[1] << "|" << str[2] << "|" << str[3] << endl;

  // TODO: Why first char should be that same as second char, followed by space?
  // TODO: Why first char should be that same as third char, delimited by space?
  if (((str[0] == str[1] && isspace(str[2])) || (str[0] == str[2] && str[1] == ' ')) && isalnum(str[0])
      && char_filter.find(str[0], 0) != string::npos)
    return str[0];

  return 0;
}

char get_atom_label(const Magick::Image &image, const Magick::ColorGray &bg, int x1, int y1, int x2, int y2,
                    double THRESHOLD, int dropx, int dropy, bool verbose)
{
  char c = 0;
  unsigned char *tmp = (unsigned char *) malloc(int((x2 - x1 + 1) * (y2 - y1 + 1)));

  for (int i = y1; i <= y2; i++)
    for (int j = x1; j <= x2; j++)
      tmp[(i - y1) * (x2 - x1 + 1) + j - x1] = (unsigned char) (255 - 255 * get_pixel(image, bg, j, i,
          THRESHOLD));

  // Here we drop down from the top of the box, middle of x coordinate and extract connected component
  int t = 1;
  int y = dropy - y1 + 1;
  int x = dropx - x1;

  while ((t != 0) && (y < int(y2 - y1 + 1)))
    {
      t = tmp[y * (x2 - x1 + 1) + x];
      y++;
    }

  if (t != 0)
    {
      free(tmp);
      return 0;
    }

  #pragma omp critical
  {
    y--;

    tmp[y * (x2 - x1 + 1) + x] = 2;

    list<int> cx;
    list<int> cy;

    cx.push_back(x);
    cy.push_back(y);

    while (!cx.empty())
      {
        x = cx.front();
        y = cy.front();
        cx.pop_front();
        cy.pop_front();
        tmp[y * (x2 - x1 + 1) + x] = 1;

        // this goes around 3x3 square touching the chosen pixel
        for (int i = x - 1; i < x + 2; i++)
          for (int j = y - 1; j < y + 2; j++)
            if ((i < (x2 - x1 + 1)) && (j < (y2 - y1 + 1)) && (i >= 0) && (j >= 0) && (tmp[j* (x2 - x1 + 1) + i] == 0))
              {
                cx.push_back(i);
                cy.push_back(j);
                tmp[j * (x2 - x1 + 1) + i] = 2;
              }
      }

    for (int i = 0; i < (y2 - y1 + 1); i++)
      for (int j = 0; j < (x2 - x1 + 1); j++)
        if (tmp[i * (x2 - x1 + 1) + j] == 1)
          {
            tmp[i * (x2 - x1 + 1) + j] = 0;
          }
        else
          {
            tmp[i * (x2 - x1 + 1) + j] = 255;
          }

    job_t gocr_job;
    double f = 1.;

    //if ((y2 - y1) > MAX_FONT_HEIGHT) f = 1. * MAX_FONT_HEIGHT / (y2 - y1);

    const int width = int(f * (x2 - x1 + 1));
    const int height = int(f * (y2 - y1 + 1));

    // The list of all characters, that can be recognised as atom label:
    string char_filter = "oOcCnNHFsSBuUgMeEXYZRPp23456789AmTh";

    job_init(&gocr_job);
    job_init_image(&gocr_job);

    //gocr_job.cfg.cs = 160;
    //gocr_job.cfg.certainty = 80;
    //gocr_job.cfg.dust_size = 1;
    gocr_job.src.p.x = width;
    gocr_job.src.p.y = height;
    gocr_job.src.p.bpp = 1;
    gocr_job.src.p.p = (unsigned char *) malloc(gocr_job.src.p.x * gocr_job.src.p.y);
    gocr_job.cfg.cfilter = (char*) char_filter.c_str();

    for (int i = 0; i < width * height; i++)
      gocr_job.src.p.p[i] = 255;

    struct OCRAD_Pixmap *ocrad_pixmap = new OCRAD_Pixmap();
    unsigned char *ocrad_bitmap = (unsigned char *) malloc(width * height);

    memset(ocrad_bitmap, 0, width * height);

    ocrad_pixmap->height = height;
    ocrad_pixmap->width = width;
    ocrad_pixmap->mode = OCRAD_bitmap;
    ocrad_pixmap->data = ocrad_bitmap;

    int count = 0;
    int zeros = 0;

    // The code below initialises the "job.src.p.p" image buffer for GOCR and "opix->data" buffer ("bitmap_data") for OCRAD from "tmp" buffer:
#ifdef HAVE_CUNEIFORM_LIB
    Magick::Image cuneiform_img(Magick::Geometry(2 * (x2 - x1 + 1) + 2, y2 - y1 + 1), "white");
    cuneiform_img.monochrome();
    cuneiform_img.type(Magick::BilevelType);
#endif
    for (int i = y1; i <= y2; i++)
      {
        for (int j = x1; j <= x2; j++)
          {
            int x = int(f * (j - x1));
            int y = int(f * (i - y1));
            if ((x < width) && (y < height) && (gocr_job.src.p.p[y * width + x] == 255))
              {
                gocr_job.src.p.p[y * width + x] = tmp[(i - y1) * (x2 - x1 + 1) + j - x1];
                if (tmp[(i - y1) * (x2 - x1 + 1) + j - x1] == 0)
                  {
                    ocrad_bitmap[y * width + x] = 1;
#ifdef HAVE_CUNEIFORM_LIB
                    cuneiform_img.pixelColor(x, y, "black");
                    cuneiform_img.pixelColor(x + (x2 - x1 + 1) + 2, y, "black");
#endif
                    if (x > 0 && x < width - 1 && y > 0 && y < height - 1)
                      count++;
                  }
                else if (x > 0 && x < width - 1 && y > 0 && y < height - 1)
                  zeros++;
              }
          }
      }

    if (verbose)
      {
        cout << "Box to OCR: " << x1 << "x" << y1 << "-" << x2 << "x" << y2 << " w/h:" << x2 - x1 << "/" << y2 - y1 << endl;
        for (int i = 0; i < height; i++)
          {
            for (int j = 0; j < width; j++)
              cout << gocr_job.src.p.p[i * width + j] / 255;
            cout << endl;
          }
      }

    if (count <= MIN_CHAR_POINTS || zeros <= MIN_CHAR_POINTS)
      goto FINALIZE;

    c = osra_gocr_ocr(gocr_job);

    if (verbose)
      cout << "GOCR: c=" << c << endl;

    //c = 0; // Switch off GOCR recognition

    // Character recognition succeeded for GOCR:
    if (c != 0)
      goto FINALIZE;

    // Character recognition failed for GOCR and we try OCRAD:
    c = osra_ocrad_ocr(ocrad_pixmap, char_filter);

    if (verbose)
      cout << "OCRAD: c=" << c << endl;

    //c = 0;  // Switch off OCRAD recognition

    // Character recognition succeeded for OCRAD:
    if (c != 0)
      goto FINALIZE;

#ifdef HAVE_TESSERACT_LIB
    c = osra_tesseract_ocr(gocr_job.src.p.p, x1, y1, x2, y2, char_filter);

    if (verbose)
      cout << "Tesseract: c=" << c << endl;

    //c = 0;  // Switch off Tesseract recognition

    // Character recognition succeeded for Tesseract:
    if (c != 0)
      goto FINALIZE;
#endif
#ifdef HAVE_CUNEIFORM_LIB
    // TODO: Why box width should be more than 7 for Cuneiform?
    // TODO: Can we replace this with "width"?
    if (x2 - x1 <= 7)
      goto FINALIZE;

    c = osra_cuneiform_ocr(cuneiform_img, verbose, char_filter);

    if (verbose)
      cout << "Cuneiform: c=" << c << endl;

    //c = 0; // Switch off Cuneiform recognition
#endif

FINALIZE:
    job_free_image(&gocr_job);
    OCR_JOB = NULL;
    JOB = NULL;

    free(tmp);
    delete ocrad_pixmap; // delete OCRAD Pixmap
    free(ocrad_bitmap);

    // TODO: Why there are problems with "7" with a given box size? If the problem is engine-specific, it should be moved to appropriate section
    // TODO: Can we replace this with "width" and "height"?
    if (c == '7' && (x2 - x1 <= 10 || y2 - y1 <= 20))
      c = 0;
  } // #pragma omp critical

  return c;
}

/*
bool detect_bracket(int x, int y, unsigned char *pic) {
	Control control;
	char c1 = 0;
	job_t job;
	JOB = &job;
	job_init(&job);
	job.cfg.cfilter = (char *) "([{";

	//job.cfg.cs = 160;
	//job.cfg.certainty = 80;
	//job.cfg.dust_size = 1;

	bool res = false;

	job.src.p.x = x;
	job.src.p.y = y;
	job.src.p.bpp = 1;
	job.src.p.p = pic;

	Blob *b = new Blob(0, 0, job.src.p.x, job.src.p.y);

	int count = 0;
	int zeros = 0;
	for (int i = 0; i <= y; i++)
		for (int j = 0; j <= x; j++) {
			if (pic[i * x + j] == 0) {
				b->set_bit(y, x, true);
				count++;
			} else
				zeros++;
		}

	if (count > MIN_CHAR_POINTS && zeros > MIN_CHAR_POINTS) {
		try {
			pgm2asc(&job);
		} catch (...) {
		}
		char *l;
		l = (char *) job.res.linelist.start.next->data;
		if (l != NULL)
			c1 = l[0];
		if (c1 == '(' || c1 == '[' || c1 == '{')
			res = true;
		else {
			char c2 = 0;
			b->find_holes();
			Character a(b);
			a.recognize1(control.charset, Rectangle::Rectangle(a.left(), a.top(), a.right(), a.bottom()));
			c2 = a.byte_result();
			if (c2 == '(' || c2 == '[' || c2 == '{')
				res = true;
		}
	}

	job_free(&job);
	return (res);
}
*/

const string fix_atom_name(const string &s, int n, const map<string, string> &fix,
                           const map<string, string> &superatom, bool debug)
{
  string r = s;

  if (s.length() == 1)
    r = toupper(s.at(0));
  if (s == "H" && n > 1)
    r = "N";

  map<string, string>::const_iterator it = fix.find(s);
  string mapped = " ";
  if (it != fix.end())
    {
      r = it->second;
      mapped = r;
    }

  if (debug && s != " " && s != "")
    {
      it = superatom.find(r);
      string smiles = " ";
      if (it != superatom.end())
        smiles = it->second;
      cout << s << " --> " << mapped << " --> " << smiles << endl;
    }

  return (r);
}
