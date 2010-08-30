// header: osra.h
// defines types and functions used globally

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

#include <limits.h>
#include <iostream>
#include <vector>

#include <Magick++.h>

extern "C" {
#include <potracelib.h>
}

using namespace std;

// struct: atom_s
// Contains information about perspective atom
struct atom_s
{
// doubles: x,y
// coordinates within the image clip
  double x, y;
// string: label
// atomic label
  string label;
// int: n
// counter of created OBAtom objects in <get_smiles()>
  int n;
// int: anum
// atomic number
  int anum;
// pointer: curve
// pointer to the curve found by Potrace
  const potrace_path_t *curve;
// bools: exists, corner, terminal
// atom exists, atom is at the corner (has two bonds leading to it), atom is a terminal atom
  bool exists, corner, terminal;
// int: charge
// electric charge on the atom
  int charge;
};
// typedef: atom_t
// defines atom_t type based on atom_s struct
typedef struct atom_s atom_t;

// struct: bond_s
// Contains information about perspective bond between two atoms
struct bond_s
{
  // ints: a,b, type
  // starting atom, ending atom, bond type (single/doouble/triple)
  int a, b, type;
  // pointer: curve
  // pointer to the curve found by Potrace
  const potrace_path_t *curve;
  // bools: exists, hash, wedge, up, down, Small, arom
  // bond existence and type flags
  bool exists;
  bool hash;
  bool wedge;
  bool up;
  bool down;
  bool Small;
  bool arom;
  // bool: conjoined
  // true for a double bond which is joined at one end on the image
  bool conjoined;
};
// typedef: bond_t
// defines bond_t type based on bond_s struct
typedef struct bond_s bond_t;

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

//struct: point_s
// a point of the image, used by image segmentation routines
struct point_s
{
  // int: x,y
  // coordinates of the image point
  int x, y;
};
// typedef: point_t
// defines point_t type based on point_s struct
typedef struct point_s point_t;

//struct: box_s
//encompassing box structure for image segmentation
struct box_s
{
  //int: x1,y1,x2,y2
  // coordinates of top-left and bottom-right corners
  int x1, y1, x2, y2;
  //array: c
  //vector of points in the box
  vector<point_t> c;
};
//typedef: box_t
//defines box_t type based on box_s struct
typedef struct box_s box_t;

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

//struct: dash_s
// used to identify dashed bonds in <osra.cpp::find_dashed_bonds()> and small bonds in <osra.cpp::find_small_bonds()>
struct dash_s
{
  //double: x,y
  //coordinates
  double x, y;
  //bool: free
  // is this dash available for a perspective dashed bond?
  bool free;
  //pointer: curve
  // pointer to the curve found by Potrace
  const potrace_path_t *curve;
  //int: area
  // area occupied by the dash
  int area;
};
//typedef: dash_t
//defines dash_t type based on dash_s struct
typedef struct dash_s dash_t;

//struct: fragment_s
// used by <osra.cpp::populate_fragments()> to split chemical structure into unconnected molecules.
struct fragment_s
{
  //int: x1,y1,x2,y2
  //top left and bottom right coordinates of the fragment
  int x1, y1, x2, y2;
  //array: atom
  //vector of atom indices for atoms in a molecules
  vector<int> atom;
};
//typedef: fragment_t
//defines fragment_t type based on fragment_s struct
typedef struct fragment_s fragment_t;


// Section: Functions
//
//  Function: fix_atom_name()
//
//  Corrects common OCR errors by using spelling dictionary
//
//  Parameters:
//
//      s - Original atomic label as returned by OCR engine.
//      n - The number of bonds attached to the atom.
//      fix - spelling dictionary
//      superatom - dictionary of superatom labels mapped to SMILES
//      debug - enables output of debugging information to stdout
//
//   Returns:
//
//      Corrected atomic label.
const string fix_atom_name(const string &s, int n, const map<string, string> &fix,
                           const map<string, string> &superatom, bool debug);

//  Function: get_smiles()
//
//  Converts vectors of atoms and bonds into a molecular object and returns SMILES, MOL file or other
//  molecular representation
//
//  Parameters:
//
//   atom - vector of <atom_s> atoms
//   bond - vector of <bond_s> bonds
//   n_bond - total number of bonds
//   rotors - number of rotatable bonds (returned)
//   confidence - confidence score (returned)
//   num_fragments - number of fragments (returned)
//   r56 - number of 5- and 6-member rings
//   avg - average bond length as measured from the image
//   format - format for molecular representation - i.e. SMI, SDF
//   resolution - resolution at which image is being processed, dpi
//   conf - toggles confidence score inclusion into output
//   guess - toggles inclusion of estimates resoltuion into output
//   showpage - toggles page number inclusion into output
//   page - page number
//   superatom - dictionary of superatom labels mapped to SMILES
//   showbond - toggles average bond length inclusion into output
//
//   Returns:
//
//    string containing SMILES, SDF or other representation of the molecule
const string get_smiles(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, int &rotors, double &confidence,
                        int &num_fragments, int &r56, double avg, const string &format, int resolution, bool conf, bool guess,
                        bool showpage, int page, const map<string, string> &superatom, bool showbond);

// Function: anisotropic_smoothing()
//
// Performs Greycstoration anisotropic smoothing on an image according to the specified parameters
//
// Parameters:
//
// image - image object
// width - width of image
// height - height of image
// amplitude - amplitude of smoothing
// alpha - alpha parameter for smoothing
// sigma - sigma parameter for smoothing'
//
// Returns:
//
// image object
//
// See also:
// <anisotropic_scaling()>
Magick::Image anisotropic_smoothing(const Magick::Image &image, int width, int height, const float amplitude,
                                    const float alpha, const float sigma);

// Function: anisotropic_scaling()
//
// Performs Greycstoration anisotropic scaling on an image
//
// Parameters:
//
// image - image object
// width - width of image
// height - height of image
// nw - new width
// nh - new height
//
// Returns:
//
// image object
//
// See also:
// <anisotropic_smoothing()>
Magick::Image anisotropic_scaling(const Magick::Image &image, int width, int height, int nw, int nh);

// Function: get_atom_label()
//
// OCR engine function, does single character recognition
//
// Parameters:
//
// image - image object
// bg - gray-level background color
// x1,y1,x2,y2 - coordinates of the character box
// THRESHOLD - graylevel threshold for image binarization
// dropx,dropy - coordinates of drop point from where breadth-first algorithm will search for single connected component
// which is hopefully the character we are trying to recognize
//
// Returns:
//
// recognized character or 0
char get_atom_label(const Magick::Image &image, const Magick::ColorGray &bg, int x1, int y1, int x2, int y2,
                    double THRESHOLD, int dropx, int dropy);

// Function: getPixel()
//
// Returns a binarized pixel value from a gray-level image
//
// Parameters:
//
// image - image object
// bg - gray-level background color
// x,y - coordinates of the pixel
// THRESHOLD - gray-level threshold for binarization
//
// Returns:
//
// 1 for set pixel, 0 for background
int getPixel(const Magick::Image &image, const Magick::ColorGray &bg, unsigned int x, unsigned int y, double THRESHOLD);

// Function: confidence_function()
//
// Calculates confidence estimate based on molecular counts provided by <get_smiles()>
//
// Parameters:
//
// C_Count, N_Count, O_Count, F_Count, S_Count, Cl_Count, Br_Count - number of carbon, nitrogen, oxygen, fluorine, sulfur, chlorine, and bromine atoms
// R_Count - number of recognized Markush atomic labels, such as R1, R2....
// Xx_Count - number of unrecognized atomic labels from <osra.cpp::remove_small_terminal_bonds()>
// num_rings - number of rings
// num_aromatic - number of aromatic rings
// num_fragments - number of fragments
// Num_Rings - vector of counts for number of 3,4,5,6,7-member rings
// num_double - number of double bonds
// num_triple - number of triple bonds
//
// Returns:
//
// confidence estimate
double confidence_function(int C_Count, int N_Count, int O_Count, int F_Count, int S_Count, int Cl_Count, int Br_Count,
                           int R_Count, int Xx_Count, int num_rings, int num_aromatic, int num_fragments, const vector<int> &Num_Rings,
                           int num_double, int num_triple);

//bool detect_bracket(int x, int y, unsigned char *pic);

// Function: unpaper()
//
// Performs unpaper image adjustment based on http://unpaper.berlios.de/
//
// Parameters:
//
// picture - image object
void unpaper(Magick::Image &picture);

// Section: Constants
//
// Constants: global defines
//
// OSRA_VERSION  - version of the program
// MAX_ATOMS  - maximum size of the vector holding perspective atoms
// MAX_FONT_HEIGHT - maximum font height at a resolution of 150 dpi
// MAX_FONT_WIDTH - maximum font width at a resolution of 150 dpi
// MIN_FONT_HEIGHT - minimum font height
// BG_PICK_POINTS - number of points to randomly pick to determine background color
// D_T_TOLERANCE - cosine tolerance to find parallel bonds for double-triple bond extraction
// V_DISPLACEMENT - threshold vertical displacement in pixels
// DIR_CHANGE - threshold direction change in pixels
// THRESHOLD_GLOBAL - gray-level threshold for image binarization
// THRESHOLD_LOW_RES - gray-level threshold for low resolutions (72 dpi)
// MAX_RATIO - maximum black/white fill ratio for perspective molecular structures
// MIN_ASPECT - minimum aspect ration
// MAX_ASPECT - maximum aspect ratio
// MIN_A_COUNT - minimum number of atoms
// MAX_A_COUNT - maximum number of atoms
// MIN_CHAR_POINTS - minimum number of black and white pixels in a character box
// MAX_BOND_THICKNESS - maximum bond thickness
// SMALL_PICTURE_AREA - threshold area of the image to be consider a small picture
// NUM_RESOLUTIONS - number of resolutions to try
// MAX_DASH - maximum size of a dash in a dashed bond
// CC_BOND_LENGTH - average carbon-carbon bond length
// FRAME - border around structure in a segmented image
// SEPARATOR_ASPECT - aspect ratio for a perspective separator line
// SEPARATOR_AREA - area for a perspective separator line
// MAX_DIST - maximum distance in pixels between neighboring segments in image segmentation routines
// MAX_AREA_RATIO - maximum area ratio for connected compoments in image segmentation
// SINGLE_IMAGE_DIST - default distance between connected components in a single structure image
// THRESHOLD_LEVEL - threshold level for feature matrix for image segmentation
// TEXT_LINE_SIZE - maximum atomic label size in characters
// PARTS_IN_MARGIN - take only every other pixel on a connected component margin for speed
// BORDER_COUNT - threshold number of pixels on a box border to be considered a table
// MAX_SEGMENTS - maximum number of connected compoment segments
// MAX_FRAGMENTS - maximum number of fragments
// STRUCTURE_COUNT - threshold number of structures to compute limits on average bond length
// SPELLING_TXT - spelling file for OCR corrections
// SUPERATOM_TXT - superatom file for mapping labels to SMILES
#define OSRA_VERSION "1.3.7"
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
#define MAX_A_COUNT 250
#define MIN_CHAR_POINTS 2
#define MAX_BOND_THICKNESS 10
#define SMALL_PICTURE_AREA 6000
#define NUM_RESOLUTIONS 4
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
#define PARTS_IN_MARGIN 2
#define BORDER_COUNT 300
#define MAX_SEGMENTS 10000
#define MAX_FRAGMENTS 10
#define STRUCTURE_COUNT 20
#define SPELLING_TXT "spelling.txt"
#define SUPERATOM_TXT "superatom.txt"
