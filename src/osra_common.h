/******************************************************************************
 OSRA: Optical Structure Recognition Application

 This is a U.S. Government work (2007-2011) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2011 (igorf@helix.nih.gov)

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
// Header: osra_common.h
//
// Common functionality routines
//

extern "C" {
#include <potracelib.h>
}

#include "osra.h"

#define BM_WORDSIZE ((int)sizeof(potrace_word))
#define BM_WORDBITS (8*BM_WORDSIZE)
#define BM_HIBIT (((potrace_word)1)<<(BM_WORDBITS-1))
#define bm_scanline(bm, y) ((bm)->map + (y)*(bm)->dy)
#define bm_index(bm, x, y) (&bm_scanline(bm, y)[(x)/BM_WORDBITS])
#define bm_mask(x) (BM_HIBIT >> ((x) & (BM_WORDBITS-1)))
#define bm_range(x, a) ((int)(x) >= 0 && (int)(x) < (a))
#define bm_safe(bm, x, y) (bm_range(x, (bm)->w) && bm_range(y, (bm)->h))
#define BM_USET(bm, x, y) (*bm_index(bm, x, y) |= bm_mask(x))
#define BM_UCLR(bm, x, y) (*bm_index(bm, x, y) &= ~bm_mask(x))
#define BM_UPUT(bm, x, y, b) ((b) ? BM_USET(bm, x, y) : BM_UCLR(bm, x, y))
#define BM_PUT(bm, x, y, b) (bm_safe(bm, x, y) ? BM_UPUT(bm, x, y, b) : 0)

//
// Section: Functions
//

// Function: get_pixel()
//
// Returns a binarized pixel value from a gray-level image
//
// Parameters:
//      image - image object
//      bg - gray-level background color
//      x, y - coordinates of the pixel
//      THRESHOLD - gray-level threshold for binarization
//
// Returns:
//      1 for set pixel, 0 for background
int get_pixel(const Magick::Image &image, const Magick::ColorGray &bg, unsigned int x, unsigned int y, double THRESHOLD);

// Function: trim()
//
// Remove leading and trailing whitespace
//
// Parameters:
//      s - string to trim (in/out parameter)
void trim(std::string &s);

double distance(double x1, double y1, double x2, double y2);
double bond_length(const vector<bond_t> &bond, int i, const vector<atom_t> &atom);
void delete_curve(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond,
                  const potrace_path_t * const curve);
void delete_curve_with_children(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond,
                                const potrace_path_t * const p);
bool detect_curve(vector<bond_t> &bond, int n_bond, const potrace_path_t * const curve);
bool terminal_bond(int a, int b, const vector<bond_t> &bond, int n_bond);
potrace_bitmap_t *const bm_new(int w, int h);
double angle4(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
double angle_between_bonds(const vector<bond_t> &bond, int i, int j, const vector<atom_t> &atom);
double distance_from_bond_y(double x0, double y0, double x1, double y1, double x, double y);
double distance_between_bonds(const vector<bond_t> &bond, int i, int j, const vector<atom_t> &atom);
double distance_from_bond_x_a(double x0, double y0, double x1, double y1, double x, double y);
double distance_from_bond_x_b(double x0, double y0, double x1, double y1, double x, double y);
double percentile75(const vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom);
int count_pages(const string &input);
int count_atoms(const vector<atom_t> &atom, int n_atom);
int count_bonds(const vector<bond_t> &bond, int n_bond, int &bond_max_type);
bool load_config_map(const string &file, map<string, string> &out);
