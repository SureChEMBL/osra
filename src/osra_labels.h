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

// Header: osra_labels.h
//
// declares functions dealing with atomic labels
//
#ifndef OSRA_LABELS_H
#define OSRA_LABELS_H

#include <string>
#include <vector>
#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX

#include <Magick++.h>

extern "C" {
#include <potracelib.h>
}


#include "osra.h"

using namespace std;

// struct: letters_s
// character found as part of atomic label
struct letters_s
{
  // double: x,y,r
  // coordinates of the center and radius of the circumscribing circle
  double x, y, r;
  // char: a
  // character
  char a;
  // bool: free
  // whether or not it was already assign to an existing atomic label
  bool free;
};
// typedef: letters_t
// defines letters_t type based on letters_s struct
typedef struct letters_s letters_t;

// struct: label_s
// atomic label
struct label_s
{
  // doubles: x1,y1, x2, y2, r1, r2
  // central coordinates and circumradii for the first and last characters
  double x1, y1, r1, x2, y2, r2;
  // string: a
  // atomic label string
  string a;
  // array: n
  // vector of character indices comprising the atomic label
  vector<int> n;
};
// typedef: label_t
// defines label_t type based on label_s struct
typedef struct label_s label_t;

//struct: lbond_s
//pairs of characters used for constucting atomic labels in <osra.cpp::assemble_labels()>
struct lbond_s
{
  //int: a,b
  // indices of first and second character in a pair
  int a, b;
  //double: x
  // x-coordinate of the first character
  double x;
  //bool: exists
  //pair of characters is available
  bool exists;
};
//typedef: lbond_t
//defines lbond_t type based on lbond_s struct
typedef struct lbond_s lbond_t;

int assemble_labels(vector<letters_t> &letters, int n_letters, vector<label_t> &label);
int find_chars(const potrace_path_t * p, const Image &orig, vector<letters_t> &letters, vector<atom_t> &atom, vector<
               bond_t> &bond, int n_atom, int n_bond, int height, int width, ColorGray &bgColor, double THRESHOLD,
               int max_font_width, int max_font_height, int &real_font_width, int &real_font_height, bool verbose);
int find_plus_minus(const potrace_path_t *p, vector<letters_t> &letters, vector<atom_t> &atom, vector<bond_t> &bond,
                    int n_atom, int n_bond, int height, int width, int max_font_height, int max_font_width, int n_letters);
int clean_unrecognized_characters(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, int real_font_height,
                                  int real_font_width, unsigned int size, vector<letters_t> &letters, int n_letters);
void remove_small_terminal_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double avg);
int remove_small_bonds(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, vector<letters_t> &letters,
                       int n_letters, int max_font_height, int min_font_height, double avg);
int find_fused_chars(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, vector<letters_t> &letters, int n_letters,
                     int max_font_height, int max_font_width, char dummy, const Image &orig, const ColorGray &bgColor,
                     double THRESHOLD, unsigned int size, bool verbose);
#endif
