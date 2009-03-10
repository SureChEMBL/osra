/*********************************************************************
  OSRA: Optical Structure Recognition

  This is a U.S. Government work (year) and is therefore not subject to copyright.  
  However, portions of this work were obtained from a GPL or GPL-compatiple source.   
  Created by Igor Filippov, 2007-2008 (igorf@helix.nih.gov)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
  USA

*********************************************************************/

#include <Magick++.h>
#include <iostream>
#include <vector>


using namespace std;
extern "C" {
#include "potracelib.h"
}

#include <limits.h>

struct atom_s {
  double x,y;
  string label;
  int n;
  potrace_path_t *curve;
  bool exists,corner;
  int charge;
};
typedef struct atom_s atom_t;

struct bond_s {
  int a,b,type;
  potrace_path_t *curve;
  bool exists;
  bool hash;
  bool wedge;
  bool up;
  bool down;
  bool Small;
  bool arom;
  bool conjoined;
};
typedef struct bond_s bond_t;

struct letters_s {
  double x,y,r;
  char a;
  bool free;
};
typedef struct letters_s letters_t;

struct label_s {
  double x1,y1,r1,x2,y2,r2;
  string a;
};
typedef struct label_s label_t;

struct point_s {
  int x,y;
};
typedef struct point_s point_t;

struct box_s {
  int x1,y1,x2,y2;
  vector<point_t> c;
};
typedef struct box_s box_t;

struct lbond_s{
  int a,b;
  double x;
  bool exists;
};
typedef struct lbond_s lbond_t;

struct dash_s {
  double x,y;
  bool free;
  potrace_path_t *curve;
  int area;
};
typedef struct dash_s dash_t;

struct fragment_s {
  int x1,y1,x2,y2;
  vector<int> atom;
};
typedef struct fragment_s fragment_t;


string fix_atom_name(string s,int n);
string get_smiles(vector<atom_t> &atom,vector<bond_t> &bond, int n_bond, int &rotors, double &confidence, int &num_fragments, int &r56,double avg,string format,int resolution,bool conf, bool guess);
Magick::Image anisotropic_smoothing(Magick::Image image,int width,int height, const float amplitude,const float alpha, const float sigma);
Magick::Image anisotropic_scaling(Magick::Image image,int width,int height, int nw, int nh);
char get_atom_label(Magick::Image image, Magick::ColorGray bg, int x1, int y1, int x2, int y2, double THRESHOLD, int dropx, int dropy);
int getPixel(Magick::Image image, Magick::ColorGray bg,unsigned int x, unsigned int y, double THRESHOLD);
double confidence_function(int C_Count,int N_Count,int O_Count,int F_Count,
			   int S_Count,int Cl_Count,int Br_Count,
			   int num_rings,int num_aromatic,
			   int num_fragments,vector<int> *Num_Rings);


#define OSRA_VERSION "1.2.1"
#define MAX_ATOMS 10000
#define MAX_FONT_HEIGHT 22
#define MAX_FONT_WIDTH 21
#define MIN_FONT_HEIGHT 5
#define BG_PICK_POINTS 100
#define D_T_TOLERANCE 0.95
#define V_DISPLACEMENT 3
#define DIR_CHANGE 2
#define THRESHOLD_GLOBAL 0.4
#define THRESHOLD_LOW_RES 0.2
#define MAX_RATIO 0.2
#define MIN_ASPECT 0.1
#define MAX_ASPECT 10.
#define MIN_A_COUNT 5
#define MAX_A_COUNT 200
#define MIN_CHAR_POINTS 2
#define WHITE_SPACE_FRACTION 0.3 
#define MAX_BOND_THICKNESS 10
#define SMALL_PICTURE_AREA 6000
#define NUM_RESOLUTIONS 3
#define MAX_DASH 14
#define CC_BOND_LENGTH 1.5120
#define FRAME 5
#define SEPARATOR_ASPECT 100
#define SEPARATOR_AREA 300
#define MAX_DIST 50
#define MAX_AREA_RATIO 50
#define SINGLE_IMAGE_DIST 100
#define THRESHOLD_LEVEL 4
#define TEXT_LINE_SIZE 8


