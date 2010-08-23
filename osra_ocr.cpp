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

extern "C" {
#include "pgm2asc.h"
}

#include <ocradlib.h>

#ifdef TESSERACT_ENABLE
#include <tesseract/baseapi.h>
#endif

#include "osra.h"

#ifdef CUNEIFORM_ENABLE
#include"ctiimage.h"
#include"cttypes.h"
#include"puma.h"
#include "lang_def.h"
#include "mpuma.h"
#include "compat_defs.h"
#endif

job_t *JOB;

/**
 * THRESHOLD is the graylevel binarization threshold.
 * dropx and dropy are the coordinates for the starting point from where the connected component (the image of the character) will be searched for.
 * Very often there is no bounding rectangle that would exclude all extraneous pieces without cutting into the character itself.
 * Those pieces confuse the OCR libraries quite a bit, so it's better to extract connected components (all characters that OSRA needs to resolve
 * are luckily single connected components) and leave the extra bits out.
 */

/**
 * The Tesseract code is supposed to be called only if both GOCR and OCRAD did't detect any alphanumeric character.
 * It is commented out because Tesseract seems to get a lot of false positives.
 */
char get_atom_label(const Magick::Image &image, const Magick::ColorGray &bg, int x1, int y1, int x2, int y2,
		double THRESHOLD, int dropx, int dropy) {
	char c = 0;
#pragma omp critical
	{
		unsigned char *tmp;
		job_t job;
		double f = 1.;

		job_init(&job);
		job.cfg.cfilter = (char *) "oOcCnNHFsSBuUgMeEXYZRPp23456789";

		//job.cfg.cs = 160;
		//job.cfg.certainty = 80;
		//job.cfg.dust_size = 1;
		//if ((y2 - y1) > MAX_FONT_HEIGHT) f = 1. * MAX_FONT_HEIGHT / (y2 - y1);


		job.src.p.x = int(f * (x2 - x1 + 1));
		job.src.p.y = int(f * (y2 - y1 + 1));
		job.src.p.bpp = 1;
		job.src.p.p = (unsigned char *) malloc(job.src.p.x * job.src.p.y);

		int height = job.src.p.y;
		int width = job.src.p.x;
		struct OCRAD_Pixmap *opix = new OCRAD_Pixmap();
		unsigned char *bitmap_data = (unsigned char *) malloc(width * height);
		memset(bitmap_data, 0, width * height);
		opix->height = height;
		opix->width = width;
		opix->mode = OCRAD_bitmap;
		opix->data = bitmap_data;

		// The code below initialises the "job.src.p.p" image buffer for GOCR and "opix" buffer for OCRAD.

		tmp = (unsigned char *) malloc(int((x2 - x1 + 1) * (y2 - y1 + 1)));

		for (int i = 0; i < job.src.p.x * job.src.p.y; i++)
			job.src.p.p[i] = 255;

		for (int i = y1; i <= y2; i++)
			for (int j = x1; j <= x2; j++)
				tmp[(i - y1) * (x2 - x1 + 1) + j - x1] = (unsigned char) (255 - 255 * getPixel(image, bg, j, i,
						THRESHOLD));

		int t = 1;
		int y = dropy - y1 + 1;
		int x = dropx - x1;

		while ((t != 0) && (y < int(y2 - y1 + 1))) {
			t = tmp[y * (x2 - x1 + 1) + x];
			y++;
		}
		y--;
		if (t == 0) {
			tmp[y * (x2 - x1 + 1) + x] = 2;
			list<int> cx;
			list<int> cy;
			cx.push_back(x);
			cy.push_back(y);
			while (!cx.empty()) {
				x = cx.front();
				y = cy.front();
				cx.pop_front();
				cy.pop_front();
				tmp[y * (x2 - x1 + 1) + x] = 1;
				for (int i = x - 1; i < x + 2; i++)
					for (int j = y - 1; j < y + 2; j++)
						if ((i < (x2 - x1 + 1)) && (j < (y2 - y1 + 1)) && (i >= 0) && (j >= 0) && (tmp[j
								* (x2 - x1 + 1) + i] == 0)) {
							cx.push_back(i);
							cy.push_back(j);
							tmp[j * (x2 - x1 + 1) + i] = 2;
						}
			}

			for (int i = 0; i < (y2 - y1 + 1); i++)
				for (int j = 0; j < (x2 - x1 + 1); j++)
					if (tmp[i * (x2 - x1 + 1) + j] == 1) {
						tmp[i * (x2 - x1 + 1) + j] = 0;
					} else {
						tmp[i * (x2 - x1 + 1) + j] = 255;
					}

			int count = 0;
			int zeros = 0;

#ifdef CUNEIFORM_ENABLE	
			Magick::Image bmp(Magick::Geometry(2*(x2-x1+1)+2,y2-y1+1),"white");
			bmp.monochrome(); 
			bmp.type(Magick::BilevelType);
#endif
			for (int i = y1; i <= y2; i++) {
				for (int j = x1; j <= x2; j++) {
					int x = int(f * (j - x1));
					int y = int(f * (i - y1));
					if ((x < job.src.p.x) && (y < job.src.p.y) && (job.src.p.p[y * job.src.p.x + x] == 255)) {
						job.src.p.p[y * job.src.p.x + x] = tmp[(i - y1) * (x2 - x1 + 1) + j - x1];
						if (tmp[(i - y1) * (x2 - x1 + 1) + j - x1] == 0) {
						        bitmap_data[y * job.src.p.x + x] = 1;
#ifdef CUNEIFORM_ENABLE	
							bmp.pixelColor(x, y, "black");
							bmp.pixelColor(x+(x2-x1+1)+2, y, "black");
#endif
							if (x > 0 && x < job.src.p.x - 1 && y > 0 && y < job.src.p.y - 1)
								count++;
						} else if (x > 0 && x < job.src.p.x - 1 && y > 0 && y < job.src.p.y - 1)
							zeros++;
					}

				}
			}

#ifdef CUNEIFORM_ENABLE	
			Magick::Blob blob;
			bmp.write(&blob, "DIB");
			size_t data_size = blob.length();
			char *dib = new char[data_size];
			memcpy(dib, blob.data(), data_size);
#endif

			/*
			cout << x2 - x1 << " " << y2 - y1 << endl;
			cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
			for (int i = 0; i < job.src.p.y; i++) {
				for (int j = 0; j < job.src.p.x; j++)
					cout << job.src.p.p[i * job.src.p.x + j] / 255;
				cout << endl;
			}
			*/

			string patern = job.cfg.cfilter;

			if (count > MIN_CHAR_POINTS && zeros > MIN_CHAR_POINTS) {
			        char c1 = 0;
				JOB = &job;
				try {
					pgm2asc(&job);
				} catch (...) {
				}
				char *l;
				l = (char *) job.res.linelist.start.next->data;
				if (l != NULL && strlen(l)==1)
				  c1 = l[0];
				//cout << "c1=" << c1 << endl;
				//c1='_';
				if (isalnum(c1)) // Character recognition succeeded for GOCR:
				  c = c1;
				else {
					// Character recognition failed for GOCR and we try OCRAD:
					char c2 = 0;
					string line = "_";
					OCRAD_Descriptor * const ocrdes = OCRAD_open();

					if (ocrdes && OCRAD_get_errno(ocrdes) == OCRAD_ok && OCRAD_set_image(ocrdes, opix, 0) == 0
							&& (job.src.p.y >= 10 || OCRAD_scale(ocrdes, 2) == 0) && OCRAD_recognize(ocrdes, 0) == 0) {
						c2 = OCRAD_result_first_character(ocrdes);
						if (OCRAD_result_blocks(ocrdes) >= 1 && OCRAD_result_lines(ocrdes, 0) && OCRAD_result_line(
								ocrdes, 0, 0) != 0)
							line = OCRAD_result_line(ocrdes, 0, 0);
					}
					//cout << "c2=" << c2 << endl;

					if (line.length() > 2)
						c2 = '_';
					OCRAD_close(ocrdes);
					if (patern.find(c2, 0) == string::npos)
						c2 = '_';
					//c2='_';
					if (isalnum(c2))
						c = c2;
#ifdef TESSERACT_ENABLE
					else {
					  char c3 = 0;
					  TessBaseAPI::InitWithLanguage(NULL, NULL, "eng", NULL, false, 0, NULL);
					  char* text = TessBaseAPI::TesseractRect(job.src.p.p, 1, x2 - x1 + 1, 0, 0, x2 - x1 + 1, y2 - y1 + 1);
					  TessBaseAPI::End();
					  if (text != NULL && strlen(text)==1)
					    c3 = text[0];
					  //cout<<"c3="<<c3<<endl;
					  if (patern.find(c3, 0) == string::npos)
					    c3 = '_';
					  if (isalnum(c3))
					    c = c3;
#endif
#ifdef CUNEIFORM_ENABLE
					  else
					    {
					      char c4=0;
					      char str[256];

					      if (x2-x1>5)
						{
						  PUMA_XOpen(dib, NULL);
						  PUMA_XFinalRecognition();
						  PUMA_SaveToMemory(NULL, PUMA_TOTEXT, PUMA_CODE_ASCII, str, 256);
						  PUMA_XClose();
						  if ((str[0]==str[1] && isspace(str[3])) || (str[0]==str[2] && str[1]==' '))
						    c4=str[0];
						  //cout<<c4<<endl;
						  if (patern.find(c4, 0) == string::npos)
						    c4 = '_';
						}
					      if (isalnum(c4))
						c = c4;
					    }
#endif
#ifdef TESSERACT_ENABLE
					}
#endif
				}
			

			}
			//cout << c << endl; // << "==========================" << endl;
#ifdef CUNEIFORM_ENABLE
			delete []dib;
#endif
		}
		job_free(&job);
		JOB = NULL;
		free(tmp);
		delete opix;
		free(bitmap_data);
		if (c == '7' && (x2 - x1 <= 10 || y2 - y1 <= 20))
			c = 0;
	} //#pragma omp critical

	if (isalnum(c)) {
		return (c);
	} else {
		return (0);
	}
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
		const map<string, string> &superatom, bool debug) {
	string r = s;

	if (s.length() == 1)
		r = toupper(s.at(0));
	if (s == "H" && n > 1)
		r = "N";

	map<string, string>::const_iterator it = fix.find(s);
	string mapped = " ";
	if (it != fix.end()) {
		r = it->second;
		mapped = r;
	}

	if (debug && s != " " && s != "") {
		it = superatom.find(r);
		string smiles = " ";
		if (it != superatom.end())
			smiles = it->second;
		cout << s << " --> " << mapped << " --> " << smiles << endl;
	}

	return (r);
}
